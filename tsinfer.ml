(*
  This program makes the inference of the selection coefficient, and population size, based on the Guassian approximation to the Wright-Fisher process.

  ===
  TODO:
  - make a better guess for alpha by using a good (not rough) guess for s. Make sure that the guess for alpha works
  - update the normal approximation for the beta function

  v.12:
  - Reduced the minimum alpha to 0.1

  v.11:
  - corrected the Gaussian approximation!!! Replaced x_0 by g(t, x_0) in the calculation of M(t) and sigma^2(t)


  v.10:
  - fixing the LLH calculation in case frequencies are not known exactly:
     . fixed the calculation of the covariance matrix 
     . fixed the calculation of the normalization constant for the normal approximation of the beta distribution (previously it took the 0-th time point into account whereas it should only start from the 1-st).

  v. 9:
  - corrected the gaussian approximation!!!

  v. 8:
  - do the computation of integrals of Gegenbauers using Laplace's method (approximating the binomials with gaussians, as in the approximate calculation)!

  v. 7:
  - replaced the explicit integration of I(x_0) over x_0 by the Laplace's method
  - added an option of evaluating only the neutral likelihood
  - disabled Monte-Carlo integration
  - removed alpha from sigmavec
  - dealing with the issue of scaling of the covariance matrix. When alpha is big and when there are many time points, the inversion of the covariance matrix is unstable because all entries (including diagonal) are close to 0.


  v. 6:
  - for s = 0, when alpha drops below 1000 switch to Kimura's experssion;
  - in the calculation of Kimura's likelihood, removed the possibility of absorbtion
  - optimized program structure


  v. 5:
  - finding a rough initial guess for s by computing the mean of d(nu)/dt/nu/(1-nu)
  - finding an initial guess for alpha by computing the inverse of the mean of (nu(t) - g(t,s))^2 / g(t,s) / (1 - g(t,s))


  v. 4:
  - using the Gaussian approximation for the binomials
  - adjusting the integration limits based on the Gaussian approximation 
  - make inference of alpha with s = 0 (to test for significance)
  - added the option of computing the likelihood function for alpha = infinity

  v. 3:
  - using the taylor expansion approximation
  - using Gsl_vector_flat


  v. 2:
  - shift time so that the initial time point is 0
  - exclude x_0 from the estimation (see v. 8 of the inf_pop_gen_time ms)
  - fixed the lack of alpha in the mean of the noise process


  v. 1:
  (based on v.14 of the temp_ml.ml file)

*)

open Bigarray
open Printf


let def_prec = 1e-6
and def_relerr = 0.001
and def_init_s = 0.
and def_init_alpha = 1000.
and def_step_s = 0.005
and def_step_alpha = 50.
and def_mins = -.2.
and def_maxs = 2.
and def_maxalpha = 1e8
and def_minalpha = 0.1
and def_kalpha = 10. (* this is the alpha below witch kimura's expression should be used *)
and def_maxiter = 200
and def_maxterms = 2000 (* maximum number of terms in the infinite sum of Gegenbauer polynomials *)
and def_func_eval = 30000
and def_stds = 3.
;;



(* === EXCEPTIONS ===== *)

(*
  code  = Small if alpha < minalpha 
        = Large if alpha > maxalpha
        = Iter if number of iterations exceeded
*)
type failcode = Success | Small | Large | Iter
;;

exception Minimizer_failed of failcode * (float, Bigarray.float64_elt, Bigarray.c_layout) Bigarray.Array1.t * float;;

exception Bad_input_line of string ;;

exception Bad_file of string ;;


(* ================================================================= *)
(* ======== matrix operations ============  *)
(* ================================================================= *)

type tridiag = { diag : float array ; above : float array ; below : float array }



let vec2mat_diag diag = 
  let n = Array.length diag in
  let mat = Array.make_matrix n n 0. in
    for i = 0 to n-1 do
      mat.(i).(i) <- diag.(i);
    done;
    Gsl.Matrix_flat.of_arrays mat


let vec2mat_tridiag m = 
  let n = Array.length m.diag in
  let _ = if Array.length m.above <> pred n || Array.length m.below <> pred n then failwith "vec2mat failed!" in
  let mat = Array.make_matrix n n 0. in
    for i = 0 to n-2 do
      mat.(i).(i) <- m.diag.(i);
      mat.(i).(succ i) <- m.above.(i);
      mat.(succ i).(i) <- m.below.(i);
    done;
    mat.(pred n).(pred n) <- m.diag.(pred n);
    Gsl.Matrix_flat.of_arrays mat



(* inversion of the tri-diagonal matrix A in the equaton Ax = f. A has the form
c1 b1 0  ... 0
a2 c2 b2 ... 0
0  a3 c3 ... 0
0  0  0  ... cn

This is actually not necessary because GSL has the Gsl_linalg.solve_tridiag routine
*)
let progonka a b c f =
  let n = Array.length c in
  let _ = if Array.length a <> n || Array.length b <> n || Array.length f <> n then failwith "progonka failed!" in
  let alpha = Array.make n nan
  and beta = Array.make n nan
  and x = Array.make n nan in
    alpha.(1) <- -. b.(0) /. c.(0);
    beta.(1) <- f.(0) /. c.(0);
    for i = 1 to n-2 do
      let d = a.(i) *. alpha.(i) +. c.(i) in
	alpha.(i+1) <- -. b.(i) /. d ;
	beta.(i+1)  <- ( f.(i) -. a.(i) *. beta.(i)  ) /. d
    done;
    x.(n-1) <- ( f.(n-1) -. a.(n-1) *. beta.(n-1) ) /. ( a.(n-1) *. alpha.(n-1) +. c.(n-1) );
    for i = n-2 downto 0 do
      x.(i) <- alpha.(i+1) *. x.(i+1) +. beta.(i+1)
    done;
    x
    



let invert_tridiag m =
  let n = Array.length m.diag in
  let _ = if Array.length m.above <> pred n || Array.length m.below <> pred n then failwith "invert_tridiag failed!"
  and diag = Gsl.Vectmat.vec_convert ~protect:true (`A m.diag)
  and abovediag = Gsl.Vectmat.vec_convert ~protect:true (`A m.above)
  and belowdiag = Gsl.Vectmat.vec_convert ~protect:true (`A m.below)
  and invm = Array.make_matrix n n nan in
    for i = 0 to n - 1 do
      let barr = Array.make n 0.
      and x = Gsl.Vectmat.vec_convert ~protect:false (`A (Array.make n nan)) in
	barr.(i) <- 1.;
	Gsl.Linalg.solve_tridiag ~diag ~abovediag ~belowdiag ~b:(Gsl.Vectmat.vec_convert ~protect:false (`A barr)) ~x;
	invm.(i) <- Gsl.Vectmat.to_array x
    done;
    Gsl.Matrix_flat.of_arrays invm
    




let det_tridiag m = 
  let n = Array.length m.diag in
  let _ = if Array.length m.above <> pred n || Array.length m.below <> pred n then failwith "invert_tridiag failed!" in
  let rec det_tridiag_r d = function
      -1 -> d
    | i when i = n - 1 -> det_tridiag_r m.diag.(i) (pred i)
    | i -> let x = m.diag.(i) *. d -. m.above.(i) *. m.below.(i) in det_tridiag_r x (pred i)
  in
    det_tridiag_r 0. (pred n)




(* apply function f to all elements of a vector *)
let map_vec_flat f v =
  Gsl.Vector_flat.of_array (Array.map f (Gsl.Vector_flat.to_array v))
;;



(* matrix - vector multiplication *)
let mult_mv_flat ~m ~v =
  let n1, n2 = Gsl.Matrix_flat.dims m in
  let y = Gsl.Vector_flat.create n1 in
  let _ = Gsl.Blas_flat.gemv Gsl.Blas_flat.NoTrans ~alpha:1. ~a:m ~x:v ~beta:0. ~y in
    y


let mult_mm_flat ~m1 ~m2 =
  let n11, n12 = Gsl.Matrix_flat.dims m1
  and n21, n22 = Gsl.Matrix_flat.dims m1 in
  let _ = if n12 <> n21 then failwith "mult_mm_flat : matrix dimensions do not match!" in
  let m = Gsl.Matrix_flat.create n11 n22 in
  let _ = Gsl.Blas_flat.gemm ~ta:Gsl.Blas_flat.NoTrans ~tb:Gsl.Blas_flat.NoTrans ~alpha:1. ~a:m1 ~b:m2 ~beta:0. ~c:m in
    m



let invert_flat m = Gsl.Vectmat.mat_flat (Gsl.Linalg.invert_LU ~protect:true (`MF m))






(* ============ end matrix operations ========================= *)





(* =============================================== *)
(* ======== Minimizer, root finder, MC integrator ============ *)
(* =============================================== *)


let print_state n minim iter =
  let x = Gsl.Vector.create n in
  let f = Gsl.Multimin.NoDeriv.minimum ~x minim in
  let ssval = Gsl.Multimin.NoDeriv.size minim in
    print_endline (Printf.sprintf "iter %d : s = %.5f, alpha = %.3e, f() = %.5g, ssize=%.2f" iter x.{0} x.{1} f ssval)
;;




(* Looking for the following situation :
             /
  \f(a)    _/
   \      /
    \    /
     \__/
      f(b)
  ________________
  a    b           
*)

let find_min_interval ~f ~min ~max ~delta0 ~mindelta ~deltafact ?(verbose=false) =
  let rec ff x (a, b) (fa, fb) delta = 
    let fx = f x in
      if verbose then print_endline (Printf.sprintf "%.5g\t%.5g\t\t%.5g\t%.5g\t|\t%.5g\t%.5g" a b fa fb x fx);
      if b == nan then
	if fx < fa
	then ff (x +. delta) (a, x) (fa, fx) delta
	else
	  let newdelta =
	    if delta > mindelta
	    then (delta /. 2.)
	    else let y = Gsl.Vector.of_array [|x|] in raise (Minimizer_failed (Small, y, fx))
	  in
	    ff (a +. newdelta)  (a, nan) (fa, nan) newdelta
      else if b = infinity then (a, infinity, infinity)
      else if fx > fb then (a, x, b)
      else if x >= max then ff infinity (b, x) (fb, fx) 0. else
	let newdelta = if delta < (max -. min) then delta *. deltafact else delta in
	  ff (x +. delta) (b, x) (fb, fx) newdelta
  in
    if verbose then Printf.printf "a\tb\tc\tf(a)\tf(b)\tf(c)\t|\tx\tf(x)\n";
    ff (min +. delta0) (min, nan) ((f min), nan) delta0
;;




let find_min f ~x1 ~x2 ~guess ~maxiter ~maxx ~abserr =
  (* print_endline (Printf.sprintf "START: x1 = %.5f; x2 = %.5f; guess = %.5f" x1 x2 guess); *)
  let xvec = Gsl.Vector.create 1
  and s = Gsl.Min.make Gsl.Min.GOLDENSECTION f guess x1 x2 in
  let rec proc i x = function
      true -> (x, f x);
    | false when i >= maxiter ->	
	Gsl.Vector.set xvec 0 x;
	raise (Minimizer_failed (Iter, xvec, f x))
    | false when x > maxx ->
	Gsl.Vector.set xvec 0 x;
	raise (Minimizer_failed (Large, xvec, f x))
    | _ ->
	Gsl.Min.iterate s;
	let (a, b) = Gsl.Min.interval s in
	let new_x = Gsl.Min.minimum s in
	let status = Gsl.Min.test_interval ~x_lo:a ~x_up:b ~epsabs:abserr ~epsrel:def_relerr in
	  proc (succ i) new_x status
  in
  let (r1,r2) as zz = proc 0 guess false in
   (*  printf "x1 = %.5f, x2 = %.5f, guess = %.5f, r1 = %.5f, r2 = %.5f done\n" x1 x2 guess r1 r2; *)
    zz
;;




(* computes the ML value for s given that alpha = infinity *)
let get_ml_alpha_inf mylf_s ~prec ~maxiter ~verbose ~guess =
  let (min, max, delta0) = match guess with
      0. -> (def_mins, def_maxs, 0.01)
    | s when s < 0. -> (guess *. 2., abs_float (1. *. guess), abs_float guess)
    | s  -> (-. guess *. 1. , guess *. 5., guess *. 2.)
  in
  let (s1, s2, s0) =
    (* (min, max, guess) *)
    try
      find_min_interval ~f:mylf_s ~min ~max ~delta0 ~mindelta:prec ~deltafact:1. ~verbose
    with
	Minimizer_failed (code, xbad, fbad) ->
	  let y = Gsl.Vector.of_array [|xbad.{0}; infinity|] in raise (Minimizer_failed (Small, y, fbad))
  in
  let _ = if verbose then print_endline (Printf.sprintf "Interval containing the ML value: (%.5f, %.5f), guess = %.5f\n" s1 s2 s0) in
    find_min mylf_s ~x1:s1 ~x2:s2 ~guess:s0 ~maxiter ~maxx:def_maxs ~abserr:prec
  





(* simplex method for the nll mimization *)
let find_min_2d (gf, gf_inf) ~iffreq ~start ~step_size ~maxiter ~maxalpha ~minalpha epsabs ?(verbose=false)  =
  let n = Gsl.Vector.length start in
  let x = Gsl.Vector.create n in
  let minim = Gsl.Multimin.NoDeriv.make Gsl.Multimin.NoDeriv.NM_SIMPLEX n gf  ~x:start ~step_size in
  let rec proc iter = 
    Gsl.Multimin.NoDeriv.iterate minim ;
    let status = Gsl.Multimin.NoDeriv.test_size minim epsabs in
      match status with
          true ->
(*	      Printf.printf "Minimum found at:\n" ;
	      print_state n minim iter; *)
	      let value = Gsl.Multimin.NoDeriv.minimum ~x minim in
		(x, value)
	| false when iter >= maxiter -> 
	    if verbose then print_state n minim iter;
(*	    let x = Gsl_vector.create n in *)
	    let value = Gsl.Multimin.NoDeriv.minimum ~x minim in
	      raise (Minimizer_failed (Iter, x, value))
	| false ->
	    let value = Gsl.Multimin.NoDeriv.minimum ~x minim in
	      if x.{1} > maxalpha then
		let _ = print_endline "Maximum alpha value reached. Trying to maximize the likelihood function at alpha = inf"
		and (mlsinf, fmininf) =
		  if iffreq then
		    let _ = print_endline "Since allele frequencies are treated as exact, there is no likelihood function to maximize. Will try to see if the latest value of s predicted the observed allele frequencies within the desired precision.\n"
		    and guess1 =  x.{0} and guess2 = start.{0} in
		      match gf_inf guess1, gf_inf guess2 with
			  x1, y1 when x1 = -.max_float && y1 = max_float -> (guess1, neg_infinity)
			| x1, y1 when x1 = max_float && y1 = -.max_float -> (guess2, neg_infinity)
			| x1, y1 when x1 = -.max_float && y1 = -.max_float -> failwith (Printf.sprintf "ERROR : Data is inconsistent with the assumption of exact frequencies. Tried s1 = %.10f and s2 = %.10f" guess1 guess2)
			| x1, y1 when x1 = infinity && y1 = infinity -> failwith (Printf.sprintf "Both guesses for s that I had (guess 1 = %.5g, guess 2 = %.5g) fail to predict the observed allele frequencies within desired precision. Should try fitting s with least squares method. The latest guess for s (s = %.5g) is still likely to be a good estimate." guess1 guess2 guess1)
			| x1, y1 -> failwith (Printf.sprintf "You should never see this error! x1 = %.5g, y1 = %.5g" x1 y1)
		  else
		    get_ml_alpha_inf gf_inf ~prec:epsabs ~maxiter ~verbose ~guess:x.{0}
		in		  
		  if fmininf > value then printf "The minimum value of NlogL at alpha = infinity, %.5g, is actually larger than the minimum value at alpha = %.5g, %.5g. The true minimum must be at alpha > %.5g (which is the upper limit on alpha), but will assume that it is at alpha = infinity.\n" fmininf x.{1} value x.{1};
		  x.{0} <- mlsinf;
		  x.{1} <- infinity;
		  (x, fmininf)
	      else
		if x.{1} < minalpha then raise (Minimizer_failed (Small, x, value))
	      else
		begin
		  if verbose then print_state n minim iter;
		  proc (succ iter) 
		end
  in
    proc 1
;;









let find_zero f alpha_1 alpha_2 alpha_guess max_iter abserr =
  let s = Gsl.Root.Bracket.make Gsl.Root.Bracket.BISECTION f alpha_1 alpha_2 in
  let rec proc i alpha = function
      true -> alpha
    | false when i > max_iter ->
	Printf.printf "Did not converge after %d iterations.\n" max_iter;
	alpha
    | _ ->
	Gsl.Root.Bracket.iterate s;
	let (a, b) = Gsl.Root.Bracket.interval s in
	let new_alpha = Gsl.Root.Bracket.root s in
	let status = Gsl.Root.test_interval a b abserr 0. in
	  proc (succ i) new_alpha status
  in
    proc 0 alpha_guess false
;;



(* =============================================== *)





(* =============================================== *)
(* ======== Auxilary functions  ============ *)
(* =============================================== *)


let filter_lh ~lh_method ~s ~alpha x =
(*  print_endline (Printf.sprintf "%.15e" x); *)
  match x with
      0. -> max_float
    | lh when lh = infinity -> -.max_float
    | lh when lh < 0. -> raise (Gsl.Error.Gsl_exn (Gsl.Error.EBADFUNC, (Printf.sprintf "%s : negative likelihood value at (s, alpha) = (%.5f, %.3e); L = %.5g\n" lh_method s alpha lh)) )
    | lh -> -. (log lh)
;;



let chi2cdf ~prec x =
    let ws = Gsl.Integration.make_ws 30 in
    let res = Gsl.Integration.qag (Gsl.Randist.chisq_pdf ~nu:1.) ~a:0. ~b:x ~epsabs:0. ~epsrel:def_relerr ~limit:30 Gsl.Integration.GAUSS31 ws in
      res.Gsl.Fun.res
;;


let get_second_deriv ~f ~x ~prec ~maxiter = 
  let dx = 2. *. prec in
  let fx = f x
  and fxplus = f (x +. dx)
  and fxminus = f (x -. dx) in
  let d2f = ( fxminus -. 2. *. fx  +. fxplus ) /. (dx**2.) in
(*    printf "fx = %.5g, fxplus = %.5g, fxminus = %.5g, denom =  %.5g, d2f = %.5g" fx fxplus fxminus ( fxminus -. 2. *. fx  +. fxplus ) d2f;
    read_line (); *)
    d2f



(*
let get_second_deriv ~f ~x ~prec ~maxiter = 
  let fx = f x in
  let rec lim y dx = function 
      0 -> printf "get_second_deriv : Maximum number of iterations exceeded. The value of the second derivative may not be accurate!\n"; y
    | i ->
	let y1 = ( (f (x -. dx)) -. 2. *. fx  +. (f (x +. dx)) ) /. dx**2. in
(*	  printf "dx = %.5g, y = %.5g, y1 = %.5g, diff =  %.5g" dx y y1 (abs_float (y1 -. y)); read_line (); *)
	  if abs_float (y1 -. y) <= prec then y else lim y1 (dx /. 2.) (pred i)
  in
    lim nan 0.01 maxiter
    *)



(* evaluates the log of the integral using the Laplace method *)
let eval_laplace ~f ~a ~b ~prec ~maxiter ~verbose =
  let funcneg x = -. (f x) in
    try
      let (x1, x2, guess) = find_min_interval ~f:funcneg ~min:a ~max:b ~delta0:0.1 ~mindelta:0.001 ~deltafact:2. ~verbose:false in
      let _ = if verbose then printf "Maximum of f(x) is in the interval [%.5g, %.5g], guess = %.5g\n" x1 x2 guess in 
      let (xmax, fmax) = find_min funcneg ~x1 ~x2 ~guess ~maxiter:def_maxiter ~maxx:1. ~abserr:prec in
      let _ = if verbose then printf "Maximum of the f(x) function found: x0 = %.5g\n" xmax in 
      let f2max = get_second_deriv ~f ~x:xmax ~prec ~maxiter in
(*	if verbose then printf "Maximum of the f(x) at xmax = %.5g:\tf(xmax) = %.5g, f''(xmax) = %.5g, int = %.5g\n" xmax fmax f2max (-. fmax +. ( (log (2. *. Gsl_math.pi)) -. (log (-.f2max)) ) /. 2.);  *)
	if abs_float f2max < prec then begin printf "Laplace method failed because the second derivative appears to be zero at the maximum.";  exit 1 end;
	if f2max > 0. then begin "Laplace method failed because the second derivative appears to be positive at the maximum (non-sense!)"; exit 1 end;
	-. fmax +. ( (log (2. *. Gsl.Math.pi)) -. (log (-.f2max)) ) /. 2.
    with
	Minimizer_failed (code, xmin_bad, fmin_bad) -> failwith "Could not evaluate the integral with respect to x0 using Laplace's methods: the minimum of the f(x) function was not found!"
	  



    
    
    
    


(* computes the mean of the product of two multivariate gaussians *)
let get_new_mu_sigma ~mu1 ~cov1 ~mu2 ~cov2 = 
  let cov1inv = invert_flat cov1
  and cov2inv = invert_flat cov2 in
  let cov3inv = Gsl.Matrix_flat.copy cov1inv in
  let _ = Gsl.Matrix_flat.add cov3inv cov2inv in
  let cov3 = invert_flat cov3inv
  and y1 = mult_mv_flat ~m:cov1inv ~v:mu1
  and y2 = mult_mv_flat ~m:cov2inv ~v:mu2 in
  let _ = Gsl.Vector_flat.add y1 y2 in
    (mult_mv_flat ~m:cov3 ~v:y1, cov3)
;;



(* computes the tridiagonal matrix I + \Sigma^{-1} \Sigma_b *)
let get_m_tmp cov_g_inv cov_b_diag =
  let n = Array.length cov_b_diag in
  let newdiag = Array.init n (function i -> 1. +. cov_g_inv.diag.(i) *. cov_b_diag.(i))
  and newabove = Array.init (n-1) (function i -> cov_g_inv.above.(i) *. cov_b_diag.(succ i) )
  and newbelow = Array.init (n-1) (function i -> cov_g_inv.below.(i) *. cov_b_diag.(i) ) in
    { diag = newdiag ; above = newabove ; below = newbelow }





(* computes the log of the integral of the product of two gaussians *)
let get_logz ~mu_b ~mu_g ~cov_b_diag ~cov_g ~log_det_g ~cov_g_inv =
  let l = float_of_int (Gsl.Vector_flat.length mu_b)
  and m_tmp = get_m_tmp cov_g_inv (Gsl.Vector_flat.to_array cov_b_diag) in
  let det_m_tmp = det_tridiag m_tmp in
  let log_det = (log det_m_tmp) +. log_det_g
  and m_tmp_inv = invert_tridiag m_tmp
  and mu12 = Gsl.Vector_flat.copy mu_g in
  let _ = Gsl.Vector_flat.sub mu12 mu_b
  and m_tmp2_inv = mult_mm_flat ~m1:m_tmp_inv ~m2:(vec2mat_tridiag cov_g_inv) in
  let y = mult_mv_flat ~m:m_tmp2_inv ~v:mu12 in
  let exparg = Gsl.Blas_flat.dot mu12 y in
(*    printf "exparg = %.5g, \tlog (2 pi det) = %.5g >> " exparg (l *. (log (2. *. Gsl_math.pi)) +. log_det); *)
    -. ( exparg +. l *. (log (2. *. Gsl.Math.pi)) +. log_det ) /. 2.



(*let get_z ~mu1 ~cov1 ~mu2 ~cov2 = 
  let cov12 = Gsl_matrix_flat.copy cov1
  and mu12 = Gsl_vector_flat.copy mu1 in
  let _ =
    Gsl_matrix_flat.add cov12 cov2;
    Gsl_vector_flat.sub mu12 mu2 in
  let cov12inv = invert_flat cov12
  and _ = Gsl_matrix_flat.scale cov12 (2. *. Gsl_math.pi) in
  let y = mult_mv_flat ~m:cov12inv ~v:mu12
  and coef = Gsl_linalg.det_LU ~protect:true (`MF cov12) in
  let exparg = Gsl_blas_flat.dot mu12 y in
    1. /. (sqrt coef) *. (exp ( -. exparg /. 2.) )
;;
*)


(* computes the scale constant for converting the beta distribution into the binomial distribution *)
let get_log_scale_const = Array.fold_left (fun sum n -> sum -. (log ((float_of_int n) +. 1.)) ) 0.


let multivariate_gaussian_logpdf ~mu ~cov x = 
  let x1 = Gsl.Vector_flat.copy x
  and cov1 = Gsl.Matrix_flat.copy cov in
  let _ =
    Gsl.Vector_flat.sub x1 mu in
  let cov_inv = invert_flat cov1
  and _ = Gsl.Matrix_flat.scale cov1 (2. *. Gsl.Math.pi) in
  let y = mult_mv_flat ~m:cov_inv ~v:x1
  and coef = Gsl.Linalg.det_LU ~protect:true (`MF cov1) in
  let exparg = Gsl.Blas_flat.dot x1 y in
    -. (exparg +. (log coef) ) /. 2.



(* analogous to multivariate_gaussian_logpdf but when the covariance matrix is diagonal *)
let multivariate_gaussian_diag_logpdf ~mu ~cov_d x = 
  let l = float_of_int (Gsl.Vector_flat.length mu)
  and x1 = Gsl.Vector_flat.copy x
  and x2 = Gsl.Vector_flat.copy cov_d in
  let _ =
    Gsl.Vector_flat.sub x1 mu;
    Gsl.Vector_flat.mul x2 x1
  in
  let exparg = Gsl.Blas_flat.dot x1 x2 in
  let logdet = Array.fold_left ( fun sum x -> sum +. (log x) ) 0. (Gsl.Vector_flat.to_array cov_d) in
    -. ( exparg +. l *. (log (2. *. Gsl.Math.pi)) +. logdet )  /. 2.


let gaussian_logpdf ~mu ~sigma x = 
  let exparg =  (x -. mu)**2. /. sigma**2. in
    -. ( exparg +. (log (2. *. Gsl.Math.pi)) +. 2. *. (log sigma ) )  /. 2.




let binoprod ~bvec ~nvec xvec =
  let bino i = Gsl.Randist.binomial_pdf bvec.(i) ~p:xvec.(i) ~n:nvec.(i) in
  let rec myfun_r acc = function 
	0 -> acc *. (bino 0)
      | k -> myfun_r (acc *. (bino k)) (k-1)
 in
  myfun_r 1. ((Array.length bvec) - 1)
;;











(* ======== Vector and matrix operations with bigarrays ============ *)
(* needed for Kimura's expression *)

(* point by point multiplication of two one-dimensional arrays *)
let point_mult a b = 
  let l = Array1.dim a in
    if l <> Array1.dim b then invalid_arg "point_mult : different lengths"
    else if l = 0 then Array1.of_array float64 c_layout [||] else
      begin
	let r = Array1.create float64 c_layout l in
	  for i = 0 to l - 1 do
	    r.{i} <- (a.{i} *. b.{i})
	  done;
	  r
      end
;;



(* inner product of two one-dimensional arrays *)
let vec_mult a b = 
  let l = Array1.dim a in
    if l <> Array1.dim b then invalid_arg "vec_mult : different lengths"
    else if l = 0 then 0. else
      begin
	let r = ref 0. in
	  for i = 0 to l - 1 do
	    r := !r +. (a.{i} *. b.{i})
	  done;
	  !r
      end
;;


(* matrix multiplication:
  a -- M x N array
  b -- N x 1 array
*)
let mat_mult a b = 
  let m = Array2.dim1 a
  and n = Array1.dim b in
    if m = 0 then invalid_arg "mat_mult : zero dimension "
    else if Array1.dim (Array2.slice_left a 0) <> n then invalid_arg "mat_mult : dimension don't correspond "
    else if n = 0 then Array1.of_array float64 c_layout [||]
    else
      let c = Array1.create float64 c_layout m in
	for i = 0 to m-1 do
	  c.{i} <- vec_mult (Array2.slice_left a i) b
	done;
	c
;;





(* =============================================== *)






(* =============================================== *)
(* ======== Likelihood calculation (approx) ============ *)
(* =============================================== *)



(* this function computes the mean vector for the gaussian which approximates the product of binomials *)
let get_mean_vec_b ~bvec ~nvec =
  let func b n = (b +. 1.) /. (n +. 2.) in
  let v = Aux.array_map2 func (Array.map float_of_int bvec) (Array.map float_of_int nvec) in
    Gsl.Vector_flat.of_array v
;;

(* this returns the diagonal of the covariance matrix for the gaussian which approximates the product of binomials *)
let get_cov_mat_b_diag ~bvec ~nvec =
  let len = Array.length bvec
  and bvec_f = Array.map float_of_int bvec
  and nvec_f = Array.map float_of_int nvec in
  let var i = let b = bvec_f.(i) and n = nvec_f.(i) in (b +. 1.) *. (n -. b +. 1.) /. (n +. 2.)**2. /. (n +. 3.) in
    Gsl.Vector_flat.of_array (Array.init len var)






let get_g ~x0 ~s t = match s with
    0. -> x0
  | _ -> x0 /. (x0 +. (1. -. x0) *. exp(-.s *. t))
;;


let get_M ~x0 ~s t = match s with
    0. -> 1.
  | _ -> exp(-.s *. t) /. (x0 +. (1. -. x0) *. exp(-.s *. t))**2.
;;


let get_var ~x0 ~s t = 
  let y0 = 1. -. x0 in
  let a = 2. *. x0 *. y0 *. t in match s with
      0. ->  a 
    | _ ->
	let b = x0**2. /. s *. (exp (s *. t))
	and c = y0**2. /. s *. (exp (-.s *. t))
	and d =  (y0**2. -. x0**2.) /. s
	and e = (2. +. s) *. x0 *. y0 in
	  (get_M ~x0 ~s t)**2. *. e *. (a +. b -. c +. d) 
;;

(* Kimura's expression for variance when s = 0, when t << N 
let get_var ~x0 ~s t = 
  let y0 = 1. -. x0 in
  let a = 2. *. x0 *. y0 *. t in match s with
      0. ->  x0 *. y0 *. t
    | _ ->
	let b = x0**2. /. s *. (exp (s *. t))
	and c = y0**2. /. s *. (exp (-.s *. t))
	and d =  (y0**2. -. x0**2.) /. s
	and e = (2. +. s) *. x0 *. y0 in
	  (get_M ~x0 ~s t)**2. *. e *. (a +. b -. c +. d) 
;;
*)


(* muvec is an array of functions that map xprime (state of the process at the prvious instant) into the expected value of the gaussian at the next time point *)
let get_muvec ~l ~x0 ~s tvec =
  let gvec = Array.map ( get_g ~x0 ~s ) tvec in
  let f i xprime =
    let dt = tvec.(i+1) -. tvec.(i) in
      gvec.(i+1) +. (xprime -. gvec.(i)) *. (get_M ~x0:gvec.(i) ~s dt)
  in
    Array.init l f

(*
let get_muvec ~l ~x0 ~s ~alpha tvec =
  let gvec = Array.map ( get_g ~x0 ~s ) tvec in
  let f i xprime = let dt = tvec.(i+1) -. tvec.(i) in gvec.(i+1) +. (xprime -. gvec.(i)) *. (get_M ~x0 ~s dt) /. (sqrt alpha) in
    Array.init l f
;;
*)


let get_sigmavec ~l ~x0 ~s tvec =
  let gvec = Array.map ( get_g ~x0 ~s ) tvec in
  let f i = let dt = tvec.(i+1) -. tvec.(i) in sqrt ( get_var ~x0:gvec.(i) ~s dt ) in
    Array.init l f



let get_Mvec ~l ~x0 ~s tvec =
  let gvec = Array.map ( get_g ~x0 ~s ) tvec in
  let f i = let dt = tvec.(i+1) -. tvec.(i) in get_M ~x0:gvec.(i) ~s dt in
    Array.init l f




(* obtains the g-covariance matrix, its inverse and log of its determinant *)
let get_cov_mat_g ~alpha ~sigmavec ~mvec =
  let len =  Array.length sigmavec in
  let cov_mat_inv = { diag = Array.make len nan ; above = Array.make (pred len) nan ; below = Array.make (pred len) nan }
  and log_det_g = ref (-. (float_of_int len) *. log alpha) in
  let _ =
    for i = 0 to (len-2) do
      log_det_g := !log_det_g +. 2. *. (log sigmavec.(i));
      cov_mat_inv.diag.(i) <- alpha /.sigmavec.(i)**2. +. alpha *. mvec.(i+1)**2. /. sigmavec.(i+1)**2. ;
      cov_mat_inv.above.(i) <-  -. alpha *. mvec.(i+1) /. sigmavec.(i+1)**2. ;
      cov_mat_inv.below.(i) <-  -. alpha *. mvec.(i+1) /. sigmavec.(i+1)**2. ;
    done;
    log_det_g := !log_det_g +. 2. *. (log sigmavec.(pred len));
    cov_mat_inv.diag.(pred len) <- alpha /. sigmavec.(pred len)**2. 
  in
  let cov_mat = invert_tridiag cov_mat_inv in
(*    printf "alpha = %.3g, sigma^2 = %.5g, len = %d, log det g = %.5f\n" alpha (sigmavec.(pred len)**2.) len !log_det_g;*)
    ( cov_mat , !log_det_g, cov_mat_inv )





(*
let get_cov_mat_g ~alpha ~sigmavec ~mvec =
  let len = Array.length sigmavec in
  let cov_mat_inv = Gsl_matrix_flat.create len len
  and lscl = len *. (log alpha) in
  let cov_mat =
    for i = 0 to (len-1) do
      for j = i to (min (i+1) (len-1)) do
	if i = j
	then	  
	  let x =
	    if i = len - 1
	    then 1./.sigmavec.(i)**2.
	    else 1./.sigmavec.(i)**2. +. mvec.(i+1)**2. /. alpha /. sigmavec.(i+1)**2. 
	  in
	    Gsl_matrix_flat.set cov_mat_inv i i x
	else
	  let x = -. mvec.(i+1) /. sigmavec.(i+1)**2. /. (sqrt alpha) in
	    Gsl_matrix_flat.set cov_mat_inv i j x;
	    Gsl_matrix_flat.set cov_mat_inv j i x;
      done;
      done;
      invert_tridiag cov_mat_int 
(*      Gsl_linalg.invert_LU ~protect:true (`MF cov_mat_inv) *)
  in
    Gsl_vectmat.mat_flat cov_mat
;;
*)



(*
let get_new_mean_sigma ~mu1 ~var1 ~mu2 ~var2 = 
  let mu3 = ( (mu1 /. var1) +. (mu2 /. var2)) /. (1./.var1 +. 1./.var2)
  and var3 = var1 *. var2 /. (var1 +. var2) in
    (mu3, sqrt var3)
;;
*)
    




let get_lh_freq ~prior ~l ~prec ~tvec ~nuvec ~s ~alpha =
  if alpha = infinity then
    let gvec = Array.map (get_g ~x0:nuvec.(0) ~s) tvec in
    let func g nu c = (abs_float (g -. nu) <= prec) && c in
      if Aux.fold_left2 func gvec nuvec true && (prior nuvec.(0) > 0.) then infinity else 0.
  else
    let muvec = get_muvec ~l ~x0:nuvec.(0) ~s tvec
    and sigmavec = get_sigmavec ~l ~x0:nuvec.(0) ~s tvec in
    let gausspdf i =
      let x = nuvec.(i+1) -. (muvec.(i) nuvec.(i)) in
	Gsl.Randist.gaussian_pdf x ~sigma:( sigmavec.(i) /. (sqrt alpha) )
    in
    let rec myfun_r acc = function 
	-1 -> acc *. (prior nuvec.(0)) 
      | k -> myfun_r (acc *. (gausspdf k) ) (k-1)
    in
    let lh = myfun_r 1. (l-1) in
(*      print_endline (Printf.sprintf "%.10f\t%.10f\t%.15e" s alpha lh); *)
      lh
      


let get_llh_freq ~prior ~l ~prec ~tvec ~nuvec ~s ~alpha =
  if alpha = infinity then
    let gvec = Array.map (get_g ~x0:nuvec.(0) ~s) tvec in
    let _ = Array.iteri (fun i g -> printf "t = %.1f : |g - nu| = %.4g\n" tvec.(i) (abs_float (g -. nuvec.(i)) ) ) gvec in
    let func g nu c = (abs_float (g -. nu) <= prec) && c
    in if (Aux.fold_left2 func gvec nuvec true) && (prior nuvec.(0) > 0.) then infinity else neg_infinity
  else
    let muvec = get_muvec ~l ~x0:nuvec.(0) ~s tvec
    and sigmavec = get_sigmavec ~l ~x0:nuvec.(0) ~s tvec in
    let rec myfun_r acc = function 
	-1 -> acc +. (log (prior nuvec.(0))) 
      | k ->
	  let x = nuvec.(k+1)
	  and mu = muvec.(k) nuvec.(k)
	  and sigma = sigmavec.(k) /. (sqrt alpha) in
	  let g = gaussian_logpdf ~mu ~sigma x in
(*	    print_endline (Printf.sprintf "k = %d\tmu = %.5f\talpha = %.2f\tsigma = %.5f\tx = %.10g\tg = %.5f" k mu alpha sigma x g);  *)
	    myfun_r (acc +. g ) (k-1)
    in
    let llh = myfun_r 0. (l-1) in
(*      print_endline (Printf.sprintf "%.10f\t%.10f\t%.15e" s alpha lh); *)
      llh





(*
(* this is the function in the internal multi-dimensional integral.
   Here l is the length of nvec, bvec, gvec, mvec, sigmavec *)
let get_int_func ~l ~nvec ~bvec ~gvec ~mvec ~sigmavec ~alpha ~x0 y =
  let f = function 
      0 -> gvec.(0)
    | i -> gvec.(i) +. (y.(i-1) -. gvec.(i-1)) *. mvec.(i) /. (sqrt alpha)
  in
  let muvec = Array.init l f (* mean of each gaussian *)
  and binopdf i = Gsl_randist.binomial_pdf bvec.(i) ~p:y.(i) ~n:nvec.(i) in
  let gausspdf i = Gsl_randist.gaussian_pdf (y.(i) -. muvec.(i)) ~sigma:(sigmavec.(i)) in
  let rec myfun_r acc = function 
      -1 -> acc 
    | k -> myfun_r ( acc *. (binopdf k) *. (gausspdf k) ) (k-1)
  in
    myfun_r 1. (l-1)
;;    


(* get limits of the internal multi-dimensional integral *)
let get_int_lims ~mu ~cov ~stds =
  let lo = Gsl_vector_flat.to_array mu
  and up = Gsl_vector_flat.to_array mu in
  let func i d =
    lo.(i) <- max 0. (lo.(i) -. stds *. (sqrt d) );
    up.(i) <- min 1. (up.(i) +. stds *. (sqrt d) );
  in
    Array.iteri func (Gsl_vector_flat.to_array (Gsl_matrix_flat.diagonal cov));
    (lo, up)
;;
*)


(* take the internal multidimensional integral by approximating the binomials by gaussians *)
let get_log_internal_integral ~l ~mu_b ~cov_b_diag ~tvec ~s ~alpha x0 =
  let tvec1 = Array.sub tvec 1 l in
  let muvec_g = Gsl.Vector_flat.of_array (Array.map (get_g ~x0 ~s) tvec1)
  and mvec = get_Mvec ~l ~x0 ~s tvec
  and sigmavec = get_sigmavec ~l ~x0 ~s tvec in
  let (cov_mat_g ,  log_det_g, cov_mat_g_inv) = get_cov_mat_g ~alpha ~sigmavec ~mvec in
  let logz =
    if alpha = infinity then
      (* IN FACT, here the product of binomials should be used.
	 However, it appears that the naive normal approximation to the beta distribution is rather poor, and deviates relatively (by a few percent, at least) strongly from the beta distribution. Therefore, it seems more consistent to use the normal approximation here as well. *)
      (*      binoprod ~bvec:bvec1 ~nvec:nvec1 (Gsl_vector_flat.to_array muvec_g) *)
      multivariate_gaussian_diag_logpdf ~mu:mu_b ~cov_d:cov_b_diag muvec_g 
    else (* Gaussian approximation for the Beta function *)
      get_logz ~mu_b ~mu_g:muvec_g ~cov_b_diag:cov_b_diag ~cov_g:cov_mat_g ~log_det_g ~cov_g_inv:cov_mat_g_inv
  in
    logz    




(* evaluating the log-likelihood using Laplace method *)
let get_llh ~prior ~l ~tvec ~nvec ~bvec ~prec ~maxiter ~verbose ~s ~alpha =
(*    Printf.printf "s = %.3g, alpha = %.3g : \n" s alpha ;
    flush_all (); *)
  let muvec_b = get_mean_vec_b ~nvec ~bvec
  and cov_mat_b_diag = get_cov_mat_b_diag ~nvec ~bvec in
  let mu_b_0 = Gsl.Vector_flat.get muvec_b 0
  and sigma_b_0 = Gsl.Vector_flat.get cov_mat_b_diag 0
  and mu_b = Gsl.Vector_flat.subvector muvec_b ~off:1 ~len:l
  and cov_b_diag = Gsl.Vector_flat.subvector cov_mat_b_diag ~off:1 ~len:l
  and logsc = get_log_scale_const (Array.sub nvec 1 l) in
  let func0 x0 = (gaussian_logpdf ~mu:mu_b_0 ~sigma:sigma_b_0 x0) +. (log (prior x0) ) 
  and func1 x0 = get_log_internal_integral ~l ~mu_b ~cov_b_diag ~tvec ~s ~alpha x0 in
  let func x0 = (func0 x0) +. (func1 x0)  in
  let logint = eval_laplace ~f:func ~a:prec ~b:(1. -. prec) ~prec ~maxiter ~verbose in 
(*  let logint = get_log_internal_integral ~l ~mu_b ~cov_b_diag ~tvec ~s ~alpha ((float_of_int bvec.(0))/.(float_of_int nvec.(0)) ) in *)
(*  let _ =
    Printf.printf "s = %.3g, alpha = %.3g : logsc = %.5g, logint = %.5g, f() = %.10g\n" s alpha logsc logint (logint +. logsc);
    flush_all () 
  in *)
    logint +. logsc

  






(* number of data points = l+1, their numbering goes from 0 to l *)
let likelihood_approx ~iffreq ~prior ~l ~nvec ~bvec ~tvec ~prec ~maxiter ~verbose ~s ~alpha =
  if alpha <= 0.
  then neg_infinity
  else match iffreq with
      true ->
	let nvec_f = Array.map float_of_int nvec
	and bvec_f = Array.map float_of_int bvec in
	let nuvec = Aux.array_map2 ( /. ) bvec_f nvec_f in
(*	  filter_lh ~lh_method:"approx" ~s ~alpha (get_lh_freq ~prior ~l ~prec ~tvec ~nuvec ~s ~alpha) *)
	  -.(get_llh_freq ~prior ~l ~prec ~tvec ~nuvec ~s ~alpha)
    | false ->
	-.(get_llh ~prior ~l ~tvec ~nvec ~bvec ~prec ~maxiter ~verbose ~s ~alpha)



(* =============================================== *)












(* =============================================== *)
(* ======== Likelihood calculation (Kimura) ============ *)
(* =============================================== *)


let get_nmax ~prec ~maxiter ~maxterms ~t =
  let f x = 0.25 *. (2. *. x +. 1.) *. x *. (x +. 1.) *. exp(-. x *. (x +. 1.) *. t) -. prec in
  let min_nmax = match f 1., f (float_of_int maxterms) with
      f1, _ when f1 < 0. -> 1
    | _, f2 when f2 > 0. -> failwith (Printf.sprintf "get_nmax : need to compute more than maxterms = %d terms to satisfy accuracy (%.2e). Increase maxterms or reduce accuracy. Exiting!" maxterms prec);
    | _, _ -> int_of_float (find_zero f 1. (float_of_int maxterms) (sqrt (-.log(prec) /. t) )  maxiter 0.1 )
  in
    min_nmax + 20
;;





let get_coef_1 ~abserr ~maxiter i n b iffreq =
  let b_f = float_of_int b 
  and n_f = float_of_int n in
  let nu = b_f /. n_f in
    match iffreq with
	true ->   nu *. (1. -. nu) *. Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. nu)
      | false ->
	  let mu = (b_f +. 1.) /. (n_f +. 2.)
	  and sigmasq = (b_f +. 1.) *. (n_f -. b_f +. 1.) /. (n_f +. 2.)**2. /. (n_f +. 3.) 
	  and f0 x = x *. (1. -. x) *. (Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. x)) in
	    if sqrt sigmasq < 0.5 then
	      let secderiv = get_second_deriv ~f:f0 ~x:mu ~prec:abserr ~maxiter in
		(f0 mu) +. sigmasq /. 2. *. secderiv
	    else
	      let f x = (Gsl.Randist.binomial_pdf b ~p:x ~n) *. (f0 x) in
	      let ws = Gsl.Integration.make_ws 30 in
	      let res = Gsl.Integration.qag f ~a:0. ~b:1. ~epsabs:abserr ~epsrel:def_relerr ~limit:30 Gsl.Integration.GAUSS31 ws in
		res.Gsl.Fun.res


(*
let int1 ~i ~mu ~sigma = 
  let f x =  x *. (1. -. x) *. (Gsl_sf.gegenpoly_n i 1.5 (1. -. 2. *. x)) *. (Gsl_randist.gaussian_pdf ~sigma (x -. mu)) in
  let ws = Gsl_integration.make_ws 30 in
  let res = Gsl_integration.qag f ~a:(-.5.) ~b:5. ~epsabs:abserr ~epsrel:def_relerr ~limit:30 Gsl_integration.GAUSS31 ws in
    res.Gsl_fun.res

let int2 ~i ~mu ~sigma =
  let f x = x *. (1. -. x) *. (Gsl_sf.gegenpoly_n i 1.5 (1. -. 2. *. x)) in
  let secderiv = get_second_deriv ~f ~x:mu ~prec:1e-5 ~maxiter:200 in
    (f mu) +. sigma**2. /. 2. *. secderiv
*)




let get_coef_l ~abserr ~maxiter i n b iffreq =
  let b_f = float_of_int b 
  and n_f = float_of_int n in
  let nu = b_f /. n_f in
    match iffreq with
	true -> Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. nu)
      | false ->
	  let mu = (b_f +. 1.) /. (n_f +. 2.)
	  and sigmasq = (b_f +. 1.) *. (n_f -. b_f +. 1.) /. (n_f +. 2.)**2. /. (n_f +. 3.) 
	  and f0 x = Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. x) in
	    if sqrt sigmasq < 0.5 then
	      let secderiv = get_second_deriv ~f:f0 ~x:mu ~prec:abserr ~maxiter in
		(f0 mu) +. sigmasq /. 2. *. secderiv
	    else
	      let f x = (Gsl.Randist.binomial_pdf b ~p:x ~n) *. (f0 x) in
	      let ws = Gsl.Integration.make_ws 30 in
	      let res = Gsl.Integration.qag f ~a:0. ~b:1. ~epsabs:abserr ~epsrel:def_relerr ~limit:30 Gsl.Integration.GAUSS31 ws in
		res.Gsl.Fun.res



(*
let get_coef_l_bar ~abserr i n b  = function
    true -> let nu = b /. n in nu *. ( Gsl_sf.gegenpoly_n i 1.5 (1. -. 2. *. nu) )
  | false ->
      let f x = (Gsl_randist.binomial_pdf (int_of_float b) ~p:x ~n:(int_of_float n)) *. x *. (Gsl_sf.gegenpoly_n i 1.5 (1. -. 2. *. x)) in
      let ws = Gsl_integration.make_ws 30 in
      let res = Gsl_integration.qag f ~a:0. ~b:1. ~epsabs:abserr ~epsrel:0. ~limit:30 Gsl_integration.GAUSS31 ws in
	res.Gsl_fun.res
;;
*)



let get_coef ~abserr ~maxiter i j n b iffreq =
  let b_f = float_of_int b 
  and n_f = float_of_int n in
  let nu = b_f /. n_f in
    match iffreq with
	true -> nu *. (1. -. nu) *. (Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. nu)) *. (Gsl.Sf.gegenpoly_n j 1.5 (1. -. 2. *. nu))
      | false ->
	  let mu = (b_f +. 1.) /. (n_f +. 2.)
	  and sigmasq = (b_f +. 1.) *. (n_f -. b_f +. 1.) /. (n_f +. 2.)**2. /. (n_f +. 3.) 
	  and f0 x = x *. (1. -. x) *. (Gsl.Sf.gegenpoly_n i 1.5 (1. -. 2. *. x)) *. (Gsl.Sf.gegenpoly_n j 1.5 (1. -. 2. *. x)) in
	    if sqrt sigmasq < 0.5 then
	      let secderiv = get_second_deriv ~f:f0 ~x:mu ~prec:abserr ~maxiter in
		(f0 mu) +. sigmasq /. 2. *. secderiv
	    else
	      let f x = (Gsl.Randist.binomial_pdf b ~p:x ~n) *. (f0 x) in
	      let ws = Gsl.Integration.make_ws 30 in
	      let res = Gsl.Integration.qag f ~a:0. ~b:1. ~epsabs:abserr ~epsrel:def_relerr ~limit:30 Gsl.Integration.GAUSS31 ws in
		res.Gsl.Fun.res

 



let update_g ~iffreq ~g_ref ~nmax ~l ~nvec ~bvec ~abserr ~maxiter =
  let (g1, g, gl) = !g_ref in
  let old_nmax = Array1.dim g1
  and new_g1 = Array1.create float64 c_layout nmax
  and new_gl = Array1.create float64 c_layout nmax
  and new_g = Array3.create float64 c_layout (l-2) nmax nmax in
    if old_nmax > 0
    then
      begin
	Array1.blit g1 (Array1.sub new_g1 0 old_nmax);
	Array1.blit gl (Array1.sub new_gl 0 old_nmax);

	for k = 0 to (l-3) do
	  for i = 0 to old_nmax - 1 do
	    Array1.blit (Array3.slice_left_1 g k i) (Array1.sub (Array3.slice_left_1 new_g k i) 0 old_nmax)
	  done;
	done
      end;

    for i = old_nmax to nmax-1 do
      new_g1.{i} <- get_coef_1 ~abserr ~maxiter i nvec.(0) bvec.(0) iffreq;
      new_gl.{i} <- get_coef_l ~abserr ~maxiter i nvec.(l-1) bvec.(l-1) iffreq;
    done;


    for k = 0 to (l-3) do
      for i = old_nmax to nmax-1 do
	
	for j = 0 to old_nmax-1 do
	  new_g.{k,i,j} <- get_coef ~abserr ~maxiter  i j nvec.(k+1) bvec.(k+1) iffreq;
	  new_g.{k,j,i} <- new_g.{k,i,j};
	done;
	for j = old_nmax to i do
	  new_g.{k,i,j} <- get_coef ~abserr ~maxiter i j nvec.(k+1) bvec.(k+1) iffreq;
	  new_g.{k,j,i} <- new_g.{k,i,j};
	done;
	
      done;
     done;

     g_ref := (new_g1, new_g, new_gl)



let update_g2 ~iffreq ~g_ref ~nmax ~nvec ~bvec ~abserr ~maxiter =
  let (g1, gl) = !g_ref in
  let old_nmax = Array1.dim g1
  and new_g1 = Array1.create float64 c_layout nmax
  and new_gl = Array1.create float64 c_layout nmax in
    if old_nmax > 0
    then
      begin
	Array1.blit g1 (Array1.sub new_g1 0 old_nmax);
	Array1.blit gl (Array1.sub new_gl 0 old_nmax);
      end;

    for i = old_nmax to nmax-1 do
      new_g1.{i} <- get_coef_1 ~abserr ~maxiter i nvec.(0) bvec.(0) iffreq;
      new_gl.{i} <- get_coef_l ~abserr ~maxiter i nvec.(1) bvec.(1) iffreq;
    done;

    g_ref := (new_g1, new_gl)






let get_const_coef ~tvec ~l ~nmax alpha =
  let v = Array2.create float64 c_layout (l-1) nmax in
    for k = 0 to l-2 do
      for i1 = 0 to nmax-1 do
	let i = i1 + 1 in
	let a = 4 * (2 * i + 1)
	and c = i * (i + 1)
	and b = (float_of_int (-i * (i+1) ) ) *. (tvec.(k+1) -. tvec.(k)) /. (2. *. alpha) in
	  v.{k,i1} <- (float_of_int a) /. (float_of_int c) *. ( exp b );
      done;
    done;
    v
;;


let get_const_coef2 ~tvec ~nmax alpha =
  let v = Array1.create float64 c_layout nmax in
    for i1 = 0 to nmax-1 do
      let i = i1 + 1 in
      let a = 4 * (2 * i + 1)
      and c = i * (i + 1)
      and b = (float_of_int (-i * (i+1) ) ) *. (tvec.(1) -. tvec.(0)) /. (2. *. alpha) in
	v.{i1} <- (float_of_int a) /. (float_of_int c) *. ( exp b );
    done;
    v
;;




(*
let get_const_coef_bar tvec l l_fix nmax alpha =
  let v = Array2.create float64 c_layout (l-l_fix) nmax in 
    for k = 0 to l-l_fix-1 do
      for i1 = 0 to nmax-1 do
	let i = i1 + 1 in
	let a = 4 * (2 * i + 1)
	and c = 2 * i * (i + 1)
	and b = (float_of_int (-i * (i+1) ) ) *. (tvec.(l_fix + k) -. tvec.(l_fix + k -1 )) /. (2. *. alpha) in
	  v.{k,i1} <- (-.1.) ** (float_of_int i) *. (float_of_int a) /. (float_of_int c) *. ( exp b );
      done;
    done;
    v
;;
*)



let sum_up g1 g gl v l =
  let rec sum_up_r tmp j = match j with
      0 ->
	let tmp1 = point_mult (Array2.slice_left v j) tmp in
	  vec_mult g1 tmp1
    | _ ->
	let tmp1 = point_mult (Array2.slice_left v j) tmp in
	  sum_up_r (mat_mult (Array3.slice_left_2 g (j-1)) tmp1) (j-1)
  in
    sum_up_r gl (l-2)



(*
let get_lh_nfix = sum_up
;;


(* k is the index of the last time point before potential fixation *)
let get_lh_fix nmax g1 g gl gl_bar_k v v_bar_k k  =
  let eye = Array1.create float64 c_layout nmax
  and v_tilde = Array2.create float64 c_layout (k+1) nmax in
    Array1.fill eye 1.;
    for k1 = 0 to k-1 do
      Array1.blit (Array2.slice_left v k1) (Array2.slice_left v_tilde k1);
    done;
    Array1.blit v_bar_k (Array2.slice_left v_tilde k);
    let ss = match k with
	0 -> gl_bar_k.{0}     
      | _ -> sum_up g1 g gl_bar_k v (k+1)
    in
      ss +. (sum_up g1 g eye v_tilde (k+2) )
;;
*)



let get_lh nmax (g1, g, gl) v  l = sum_up g1 g gl v l

let get_lh2 nmax (g1, gl) v = vec_mult g1 (point_mult gl v)



let likelihood_kimura ~iffreq ~nmax_ref ~g_ref ~l ~nvec ~bvec ~tvec ~mindt ~prec ~maxiter ~maxterms alpha =
  let min_nmax = get_nmax ~prec ~maxiter ~maxterms ~t:(mindt /. 2. /. alpha) in
    if min_nmax > !nmax_ref
    then begin nmax_ref := min_nmax; update_g ~iffreq ~g_ref ~nmax:(!nmax_ref) ~l ~nvec ~bvec ~abserr:prec ~maxiter; end;
    let v = get_const_coef ~tvec ~l ~nmax:(!nmax_ref) alpha in
    let lh = get_lh !nmax_ref !g_ref v l in
(*      print_endline (Printf.sprintf "alpha = %.6f ; L = %.10e ; nmax = %d ;  min_nmax = %d" alpha l !nmax_ref min_nmax);  *)
      filter_lh ~lh_method:"kimura" ~s:0. ~alpha lh


(* kimura's likelihood in case there are only 2 time points *)
let likelihood_kimura_2 ~iffreq ~nmax_ref ~g_ref ~nvec ~bvec ~tvec ~mindt ~prec ~maxiter ~maxterms alpha =
  let min_nmax = get_nmax ~prec ~maxiter ~maxterms ~t:(mindt /. 2. /. alpha) in
    if min_nmax > !nmax_ref
    then begin nmax_ref := min_nmax; update_g2 ~iffreq ~g_ref ~nmax:(!nmax_ref) ~nvec ~bvec ~abserr:prec ~maxiter; end;
    let v = get_const_coef2 ~tvec ~nmax:(!nmax_ref) alpha in
    let lh = get_lh2 !nmax_ref !g_ref v in
(*      print_endline (Printf.sprintf "alpha = %.6f ; L = %.10e ; nmax = %d ;  min_nmax = %d" alpha l !nmax_ref min_nmax);  *)
      filter_lh ~lh_method:"kimura" ~s:0. ~alpha lh
  


(* =============================================== *)


























(* =============================================== *)
(* ======== Parsing and checking data  ============ *)
(* =============================================== *)


let read_sample_file filename =
  let shift dx x = x -. dx
  and f (tvec, bvec, nvec, i, len) line =
    try
      match (Str.split (Str.regexp "\t") line), i with
	  [], _ -> (tvec, bvec, nvec, i, len)
	| l, 0 ->
	    let arr = Array.of_list (List.map float_of_string l) in
	      (arr, bvec, nvec, succ i, (Array.length arr))
	| l, 1 ->
	    let arr = Array.of_list (List.map int_of_string l) in
	      if Array.length arr = len then
		(tvec, arr, nvec, succ i, (Array.length arr))
	      else
		raise (Bad_file (Printf.sprintf "Lines of file %s must have equal lengths" filename))
	| l, 2 ->
	    let arr = Array.of_list (List.map int_of_string l) in
	      if Array.length arr = len then
		(tvec, bvec, arr, succ i, (Array.length arr))
	      else
		raise (Bad_file (Printf.sprintf "Lines of file %s must have equal lengths" filename))
	| _ -> (tvec, bvec, nvec, i, len)
    with
	Failure s -> match s with
	    "int_of_string" -> raise (Bad_file (Printf.sprintf "Second and third lines of input file %s must contain integers values" filename))
	  | "float_of_string" -> raise (Bad_file (Printf.sprintf "First line of input file %s must contain numbers" filename))
  in
  let (tvec, bvec, nvec, i, len) = List.fold_left f ([||], [||], [||], 0, 0) (Io.read_file_to_list filename) in
    (len-1, Array.map (shift tvec.(0)) tvec, bvec, nvec)
;;





let read_plot_file filename =
  let f arrl line =
    let l = Str.split (Str.regexp "\t") line in match List.length l with
	0 -> arrl
      | 2 ->
	  let arr = Array.of_list (List.map float_of_string l) in
	  let vec = Gsl.Vector.of_array arr in
	    vec::arrl
      | _ -> raise (Bad_file (Printf.sprintf "file %s has %d entries in each line instead of 2" filename (List.length l)) )
  in
  let l = List.fold_left f [] (Io.read_file_to_list filename) in
    match l with
	[] -> raise (Bad_file (Printf.sprintf "file %s is empty" filename))
      | _ -> List.rev l
;;




(* this function checks that :
   a) at least two data points
   b) there is no potential absorbtion event if iffreq = true
*)
let check_data nvec bvec iffreq = 
  if Array.length nvec > 1 then
    if iffreq then
      let func b n =
	if b = 0 || b = n then raise (Bad_file "If sampled frequencies are treated as exact, no absorbtion events are admissible")
      in
	Aux.array_iter2 func bvec nvec
    else ()
  else
    raise (Bad_file "There must be more than one data point to make estimation")
;;



(* =============================================== *)




(* =============================================== *)
(* ======== Parsing parameters  ============ *)
(* =============================================== *)
let info_string = Printf.sprintf "
TEMP_ML_SEL

  Program that makes a maximum likelihood estimate of the population size (alpha) and the selection coefficient (s) from a time-series population sample.

  (c) Sergey Kryazhimskiy 2008 -- 2014

  Reference: Feder AF, Kryazhimskiy S, Plotkin (2014). Identifying signatures of selection in time-series allele-frequency data. Genetics 196:509--522

  NOTES:
  (1) Time is assumed to be measured in generations
  (2) The default procedure of calculating the likelihood of data takes binomial sampling into account. However, this procedure was not extensively tested and may contain errors. Use instead the procedure that treats sampled frequecies as exact (-f option, see below).

  USAGE:

  ./temp_ml_sel infile [OPTIONS]

  infile -- file containing the time series data. This file should contain three lines (see examples), each tab separated, with equal number of entries. The first line contains sample times, the second line contains the number of samples of the focal allele (must be integers) and the third line containes total number of samples at each time point (must be integers).


  OPTIONS:
  -f -- treat sampled frequencies as actual frequencies
  -func_eval:XX -- number of function evalutations in the Monte Carlo integration (positve integer). Default : %.0e
  -h -- print this help
  -maxiter:XX -- maximum number of iterations (positive integer). Default : %d
  -maxterms:XX -- maximum number of terms in the series (positive integer). Default %d
  -neut -- evaluate only the likelihood of the data under neutrality
  -prec:XX -- precision (for the size of the simples). Default : %.1e
  -start:s,alpha -- starting guess for s (float) and alpha (positive float). Default : (%.5g, %.3e)
  -stds:XX -- number of standard deviations around the mean of the beta distribution for the approximation of the intergral (positive integer). Default : %.1f
  -step:s,alpha -- step size for s, and alpha (all positive floats). Default : (%.5g, %.3e)
  -plot:filename -- instead of finding the maximum likelihood value, plot the likelihood function profile at points specified in the file 'filename'. Each line of the file must consist of 2 tab-separated numbers corresponding to the coordinates (s, alpha). See examples. The 'start', and 'step' options are ignored.
  -v -- verbose mode

" (float_of_int def_func_eval) def_maxiter def_maxterms def_prec def_init_s def_init_alpha def_stds def_step_s def_step_alpha
;;


let parse_params (start_vec, step_vec, maxiter, maxterms, prec, func_eval, stds, iffreq, ifneut, mc, plot, conf, verbose) x = match (Str.split (Str.regexp ":") x) with
    ["h"] -> print_endline info_string; exit 0;
  | ["maxiter"; y] ->
      let y1 = int_of_string y in
	if y1 > 0 then maxiter := y1 else raise (Bad_input_line ("Wrong sign: "^x))
  | ["maxterms"; y] ->
      let y1 = int_of_string y in
	if y1 > 0 then maxterms := y1 else raise (Bad_input_line ("Wrong sign: "^x))
  | ["prec"; y] ->
      let y1 = float_of_string y in
	if y1 > 0. then prec := y1 else raise (Bad_input_line ("Wrong sign: "^x))
  | ["func_eval"; y] ->
      let y1 = int_of_string y in
	if y1 > 0 then func_eval := y1 else raise (Bad_input_line ("Wrong sign: "^x))
  | ["stds"; y] ->
      let y1 = float_of_string y in
	if y1 > 0. then stds := y1  else raise (Bad_input_line ("Wrong sign: "^x))
  | ["start"; y] -> begin
	try
	  let [s; alpha] = List.map float_of_string (Str.split (Str.regexp ",") y) in
	    if alpha <= 0. then raise (Bad_input_line ("Wrong sign: "^x));
	    start_vec.{0} <- s;
	    start_vec.{1} <- alpha;
	with
	    Match_failure _ -> raise (Bad_input_line ("Wrong format: "^x) )
    end;
  | ["step"; y] -> begin
	try
	  let [s; alpha] = List.map float_of_string (Str.split (Str.regexp ",") y) in
	    if s <= 0. then raise (Bad_input_line ("Wrong sign: "^x));
	    if alpha <= 0. then raise (Bad_input_line ("Wrong sign: "^x));
	    step_vec.{0} <- s;
	    step_vec.{1} <- alpha;
      with
	  Match_failure _ -> raise (Bad_input_line ("Wrong format: "^x))
    end
  | ["f"] -> if !mc then raise (Bad_input_line ("Cannot use -f and -mc options simultaneously")) else iffreq := true
  | ["neut"] -> ifneut := true
(*  | ["mc"] -> if !iffreq then raise (Bad_input_line ("Cannot use -f and -mc options simultaneously")) else mc := true *)
  | ["plot"; filename] ->
      if Sys.file_exists filename
      then plot := read_plot_file filename
      else raise (Bad_input_line ("Could not find file : "^filename))
  | ["conf"] -> conf := true
  | ["v"] -> verbose := true
  | _ -> raise (Bad_input_line ("Unknown option: "^x))
;;




let get_options options = 
  let start_vec = Gsl.Vector.of_array [|nan; -.1.|]
  and step_vec = Gsl.Vector.of_array [|nan; -.1.|]
  and maxiter = ref def_maxiter
  and maxterms = ref def_maxterms
  and prec = ref def_prec
  and func_eval = ref def_func_eval
  and stds = ref def_stds
  and iffreq = ref false
  and ifneut = ref false
  and mc = ref false
  and plot = ref []
  and conf = ref false
  and verbose = ref false in
(*    start_vec.{0} <- def_init_s;
    start_vec.{1} <- def_init_alpha;
    step_vec.{0} <- def_step_s;
    step_vec.{1} <- def_step_alpha; *)
    List.iter (parse_params (start_vec, step_vec, maxiter, maxterms, prec, func_eval, stds, iffreq, ifneut, mc, plot, conf, verbose)) options;
    (start_vec, step_vec, !maxiter, !maxterms, !prec, !func_eval, !stds, !iffreq, !ifneut, !mc, !plot, !conf, !verbose)
;;


(* =============================================== *)










let print_params start step maxiter prec func_eval stds iffreq tvec nvec bvec =
  let func_f x = Printf.sprintf "%.0f" x
  and func_int x = Printf.sprintf "%d" x in
    print_endline (Printf.sprintf "maxiter = %d
prec = %.1e
func_eval = %.0e
stds = %.1f" maxiter prec (float_of_int func_eval) stds);
    print_endline (Printf.sprintf "Initial values : (s, alpha) = (%.5f, %.3e)
Initial step size (s, alpha) = (%.5f, %.3e)" start.{0} start.{1} step.{0} step.{1});
    print_endline "\n---";
    print_string "tvec = ["; print_string (String.concat "; " (List.map func_f (Array.to_list tvec))); print_string "]\n";
    print_string "bvec = ["; print_string (String.concat "; " (List.map func_int (Array.to_list bvec))); print_string "]\n";
    print_string "nvec = ["; print_string (String.concat "; " (List.map func_int (Array.to_list nvec))); print_string "]\n";
    if iffreq then print_endline "Treating sampled frequencies as exact";
    print_endline "";
;;




let get_mindt tvec =
  let mylength = Array.length tvec
  and f (mindt, t1) t2 = if (t2 -. t1) < mindt then (t2 -. t1, t2) else (mindt, t2) in
  let (mindt, t) = Array.fold_left f (tvec.(1) -. tvec.(0), tvec.(0) ) (Array.sub tvec 1 (mylength - 1) ) in
    mindt
;;



(* auxilary funtion : computes the ML value with the s parameter constrained to s = 0 and performs the LLR test *)
let get_ml_s_0 ~maxiter ~prec ~verbose (mylf_0, mylf_k) alpha_guess =
  let delta0 = if alpha_guess = infinity then def_maxalpha/.2. else alpha_guess /. 2. in
  let _ = Printf.printf "alpha guess = %.3e\n" alpha_guess; flush_all () in
  let ((alpha1, alpha2, alpha0), myfunc) =
    try
      let interval = find_min_interval ~f:mylf_0 ~min:def_kalpha ~max:def_maxalpha ~delta0 ~mindelta:prec ~deltafact:2. ~verbose in
	(interval, mylf_0)
    with
	Minimizer_failed (code, xbad, fbad) ->
	  print_endline "WARNING! ML alpha is in the zone where the approximation breaks down. Switching to Kimura's expression."; 
	  try
	    let interval = find_min_interval ~f:mylf_k ~min:1. ~max:def_maxalpha ~delta0 ~mindelta:prec ~deltafact:2. ~verbose in
	      (interval, mylf_k)
	  with
	      Minimizer_failed (code, xbad, fbad) ->
		let y = Gsl.Vector.of_array [|0.; xbad.{0}|] in raise (Minimizer_failed (Small, y, fbad))
  in
  let _ = if verbose then print_endline (Printf.sprintf "Interval containing the ML value: (%.3e, %.3e), guess = %.3e\n" alpha1 alpha2 alpha0) in
    if alpha0 = infinity then
      (alpha0, mylf_0 infinity)
    else
      find_min myfunc ~x1:alpha1 ~x2:alpha2 ~guess:alpha0 ~maxiter ~maxx:def_maxalpha ~abserr:prec
;;




(* finds a good initial guess for s *)
(* let get_start_vec mylf_s ~prec ~maxiter ~verbose start_vec =
let get_start_vec mylf_s ~prec ~maxiter ~verbose start_vec
  let new_start_vec = Gsl_vector.copy start_vec 
  and s = start_vec.{0} in
  let _ = if s = def_init_s then
      let (mls, mlval) =
	print_endline "Finding a good initial guess for s:\n";
        get_ml_alpha_inf mylf_s ~prec ~maxiter ~verbose ~guess:def_init_s
      in
	new_start_vec.{0} <- mls
  in
    new_start_vec
;;
*)

(*  finding a rough initial guess for s by computing the mean of d(nu)/dt/nu/(1-nu)
    This follows from \dot nu = s \nu (1 - \nu)
*)
let get_rough_s_guess ~l ~tvec ~nuvec =
  let dt_vec = Aux.array_map2 ( -. ) (Array.sub tvec 1 (l-1)) (Array.sub tvec 0 (l-1))
  and dnu_vec = Aux.array_map2 ( -. ) (Array.sub nuvec 1 (l-1)) (Array.sub nuvec 0 (l-1))
  and nu_mid_vec = Aux.array_map2 ( fun a b -> (a +. b) /. 2. ) (Array.sub nuvec 1 (l-1)) (Array.sub nuvec 0 (l-1)) in 
  let dnudt_vec = Aux.array_map2 ( /. ) dnu_vec dt_vec in
  let s_est_vec = Aux.array_map2 ( fun a b -> a /. b /. (1. -. b) ) dnudt_vec nu_mid_vec in
(*  let _ =
    print_endline (String.concat "; " (List.map string_of_float (Array.to_list s_est_vec)));
  in *)
  let (sum_s, cnt) = Array.fold_left ( fun (s, c) x -> if x == nan then (s, c) else (s+.x, c+.1.) ) (0., 0.) s_est_vec in
  let s_guess = sum_s /. cnt in 
    if abs_float s_guess <= def_prec then 0. else s_guess
;;

(* finding an initial guess for alpha by computing the inverse of the mean of (nu(t) - g(t,s))^2 / g(t,s) / (1 - g(t,s))
   This follows from the fact that the variance of the noise process is
   \sigma^2 = 1/N * \nu * (1 - \nu)
*)
let get_rough_alpha_guess ~l ~tvec ~bvec ~nvec ~s =
  let nvec_f = Array.map float_of_int nvec
  and bvec_f = Array.map float_of_int bvec in
  let nuvec = Aux.array_map2 ( /. ) bvec_f nvec_f in
  let init_x = if s = 0. then (Array.fold_left ( +. ) 0. nuvec) /. (float_of_int l) else nuvec.(0) in
  let gvec = Array.map (get_g ~x0:init_x ~s) (Array.sub tvec 1 (l-1)) in
  let tmp_vec = Aux.array_map2 (fun x y -> (x -. y)**2. /. y /. (1. -. y) ) (Array.sub nuvec 1 (l-1)) gvec in
(*  let _ =
    Printf.printf "s = %.15f; l = %d; init_x = %.15f\n" s l init_x;
    print_endline (String.concat "; " (List.map string_of_float (Array.to_list nuvec)));
    print_endline (String.concat "; " (List.map string_of_float (Array.to_list gvec)));
    print_endline (String.concat "; " (List.map string_of_float (Array.to_list tmp_vec)))
  in *) 
  let (sum_inv_n, cnt) = Array.fold_left ( fun (s, c) x -> if x == nan then (s, c) else (s+.x, c+.1.) ) (0., 0.) tmp_vec in
    match sum_inv_n with
	0. -> def_maxalpha /. 2.
      | _ -> max (def_minalpha *. 2.) (cnt /. sum_inv_n)
;;
    




let get_start_vec start_vec step_vec ~tvec ~bvec ~nvec =
  let l = Array.length tvec 
  and nvec_f = Array.map float_of_int nvec
  and bvec_f = Array.map float_of_int bvec in
  let nuvec = Aux.array_map2 ( /. ) bvec_f nvec_f in
  let s_guess = get_rough_s_guess ~l ~tvec ~nuvec in
  let alpha_guess = get_rough_alpha_guess ~l ~tvec ~bvec ~nvec ~s:s_guess in
  let new_start_vec = Gsl.Vector.copy start_vec
  and new_step_vec = Gsl.Vector.copy step_vec in
(*    print_endline (Printf.sprintf "yyy%.5f\t%.5f\n" start_vec.{0} start_vec.{1}); *)
    if start_vec.{1} = -.1.
    then
      begin
	new_start_vec.{0} <- s_guess;
	new_start_vec.{1} <- alpha_guess;
      end;
    if step_vec.{1} = -.1.
    then
      begin
	new_step_vec.{0} <- max (abs_float (new_start_vec.{0} /. 10.)) def_step_s;
	new_step_vec.{1} <- max (new_start_vec.{1} /. 10.) def_step_alpha;
      end;
    (new_start_vec, new_step_vec)
;;





(* auxilary function : handles the output of the main routine *)
let output_res ~status ~prec [|mls; mlalpha; nll; mlalpha0; nll0; bads; badalpha; badnll|] =
  let llr = 2.*.(nll0 -. nll) in
  let chi2p = if llr >= 0. then 1. -. (chi2cdf ~prec llr) else nan
  and mstate = Printf.sprintf "\nminimizer status:\nS = %-.5f\nALPHA = %-.3e\nNLL = %-.5g\n" bads badalpha badnll in
  let (string1, string2) = match status with
      Success -> "SUCCESS", ""
    | Iter -> "FAILURE (Maximum number of iterations exceeded)", mstate
    | Large -> "FAILURE (Exceeded maximum alpha value and convergence at infinity failed)", mstate
    | Small -> "FAILURE (Minimum alpha value surpassed. Approximation invalid)", mstate
  in
    Printf.printf "\nSTATUS = %s
MLS = %-.5f
MLALPHA = %-.3e
MNLOGL = %-.5g
AIC2 = %-.5g
S0ALPHA = %-.3e
S0NLOGL = %-.5g
AIC1 = %-.5g
LLR = %-.5g
CHI2P = %.1e
%s" string1 mls mlalpha nll (2. *. ( 2. +. nll)) mlalpha0 nll0 (2. *. ( 1. +. nll0)) llr chi2p string2;
    exit 1
;;





 

let main () =
  try
    let (params, options) = Io.parse_input_line (Sys.argv) in
    let filename = match params with
	[] -> raise (Bad_input_line "No input file name")
      |	[x] -> x
      | _ -> raise (Bad_input_line "Unknown parameter(s)")
    in
    let (my_start_vec, my_step_vec, maxiter, maxterms, prec, calls, stds, iffreq, ifneut, mc, plot, conf, verbose) = get_options options
    and (l, tvec, bvec, nvec) = read_sample_file filename in
    let _ = check_data nvec bvec iffreq
(*    and nvec_f = Array.map float_of_int nvec
    and bvec_f = Array.map float_of_int bvec in
    let nuvec = Aux.array_map2 ( /. ) bvec_f nvec_f *)
    and prior x = 1.
    and rng = Gsl.Rng.make (Gsl.Rng.default ()) in
    let mylf ~x = likelihood_approx ~iffreq ~prior ~l ~nvec ~bvec ~tvec ~prec ~maxiter ~verbose ~s:x.{0} ~alpha:x.{1}
    and mylf_inf x = likelihood_approx ~iffreq ~prior ~l ~nvec ~bvec ~tvec ~prec ~maxiter ~verbose ~s:x ~alpha:infinity
    and mylf_0 x = likelihood_approx ~iffreq ~prior ~l ~nvec ~bvec ~tvec ~prec ~maxiter ~verbose ~s:0. ~alpha:x
    and nmax_ref = ref 0 in
    let mylf_k =
      if l > 1 then
	let g1 = Array1.create float64 c_layout 0
	and gl = Array1.create float64 c_layout 0
	and g = Array3.create float64 c_layout (l-2) !nmax_ref !nmax_ref in
	let g_ref = ref (g1, g, gl) in
	  likelihood_kimura ~iffreq ~l:(l+1) ~nvec ~bvec ~tvec ~mindt:(get_mindt tvec) ~prec ~maxiter ~maxterms ~nmax_ref ~g_ref
      else
	let g1 = Array1.create float64 c_layout 0
	and gl = Array1.create float64 c_layout 0 in
	let g_ref = ref (g1, gl) in
	  likelihood_kimura_2 ~iffreq ~nvec ~bvec ~tvec ~mindt:(get_mindt tvec) ~prec ~maxiter ~maxterms ~nmax_ref ~g_ref
	
    in
    let (start_vec, step_vec) = get_start_vec my_start_vec my_step_vec ~tvec ~bvec ~nvec in
    let _ = print_params start_vec step_vec maxiter prec calls stds iffreq tvec nvec bvec
    in
      if plot = [] then
	let status = ref Success
	and res_arr = Array.make 8 nan in
	  begin
	    try
	      if ifneut
	      then
		let alpha_guess = max def_minalpha (get_rough_alpha_guess ~l ~nvec ~bvec ~tvec ~s:0.) in
		let (mlalpha_s_0, mlval_s_0) = get_ml_s_0 ~maxiter ~prec ~verbose (mylf_0,mylf_k) alpha_guess in 
		  print_endline (Printf.sprintf "Minimum found : (s, alpha) = (0, %.3e); -log L = %.5g\n" mlalpha_s_0 mlval_s_0);
		  res_arr.(3) <- mlalpha_s_0;
		  res_arr.(4) <- mlval_s_0
	      else
		let (xmin, fmin1) = find_min_2d (mylf, mylf_inf) ~iffreq ~start:start_vec ~step_size:step_vec ~maxiter:maxiter ~maxalpha:def_maxalpha ~minalpha:def_minalpha 0.1 ~verbose in
		  res_arr.(0) <- xmin.{0};
		  res_arr.(1) <- xmin.{1};
		  res_arr.(2) <- fmin1;
		  print_endline (Printf.sprintf "Minimum found : (s, alpha) = (%.5f, %.3e); -log L = %.5g\nMaximizing log-likelihood for s = 0 (fixed)" xmin.{0} xmin.{1} fmin1);
		  if abs_float xmin.{0} <= prec then begin res_arr.(3) <- xmin.{1}; res_arr.(4) <- fmin1 end
		  else
		    let alpha_guess = max def_minalpha (get_rough_alpha_guess ~l ~nvec ~bvec ~tvec ~s:0.) in
		    let (mlalpha_s_0, mlval_s_0) = get_ml_s_0 ~maxiter ~prec ~verbose (mylf_0,mylf_k) alpha_guess in
		      print_endline (Printf.sprintf "Minimum found : (s, alpha) = (0, %.3e); -log L = %.5g\n" mlalpha_s_0 mlval_s_0);
		      res_arr.(3) <- mlalpha_s_0;
		      res_arr.(4) <- mlval_s_0
	  with
	      Minimizer_failed (code, xmin_bad, fmin_bad) ->
		res_arr.(5) <- xmin_bad.{0};
		res_arr.(6) <- xmin_bad.{1};
		res_arr.(7) <- fmin_bad;
		status := code;
(*	    | _ -> failwith "Failure!"; *)
	  end;
	  output_res ~status:(!status) ~prec res_arr
      else (* if plot = [] *)
	let func x = 
	  Printf.printf "%.5f\t%.3e\t%.5g\n" x.{0} x.{1} (mylf ~x);
	  flush_all ()
	in
	  print_endline "s\talpha\tNlogL";
	  List.iter func plot
  with
      Bad_input_line s ->
	print_endline info_string;
	print_endline (Printf.sprintf "Bad input line! %s" s); exit 1;
    | Bad_file s -> print_endline (Printf.sprintf "Bad file format : %s" s); exit 1;
    | Failure s -> print_endline (Printf.sprintf "Failure : %s" s); exit 1;
    | Gsl.Error.Gsl_exn (errno, s) -> print_endline (Printf.sprintf "Gsl error %s : %s" (Gsl.Error.string_of_errno errno) s); exit 1;
;;



Gsl.Error.init ();;
main ();;



