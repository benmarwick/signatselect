(* This module encompasses various useful auxilary functions *)

Random.init (int_of_float (Unix.time())) ;;


(* =============================================== *)
(* ======== List Functions  ==================== *)
(* =============================================== *)

let splice l o len =
  let f (l1, k) e =
    if k < o then (l1, k+1) else if List.length l1 = len then (l1, k+1) else (e::l1, k+1) in
    if List.length l < len then invalid_arg "length of the splice region is bigger than the list itself";
    if o > (List.length l) - len then invalid_arg "the splice region does not appear to be entirely within the list";
    let (new_l, _) = List.fold_left f ([], 0) l in
    List.rev new_l
    

let purge l =
  let rec f new_l old_l = match new_l, old_l with
      _, [] -> new_l
    | [], hd::tl -> f [hd] tl
    | a::tl_n, b::tl_o -> if b = a then f new_l tl_o else f (b::new_l) tl_o
  in
    List.rev (f [] l)


(* =============================================== *)
(* ======== Array Functions  ==================== *)
(* =============================================== *)

let matrix_tr mat = 
  let n1 = Array.length mat
  and n2 = Array.length mat.(0) in
  let newmat = Array.make_matrix n2 n1 mat.(0).(0) in
    for i = 0 to pred n1 do
      for j = 0 to pred n2 do
	newmat.(j).(i) <- mat.(i).(j);
      done
    done;
    newmat


let matrix_copy mat =
  let n = Array.length mat in
  let copy = Array.make n [||] in
    for i = 0 to n - 1 do
      copy.(i) <- Array.copy mat.(i)
    done;
    copy



(* fold left with index for arrays *)
let fold_lefti f x a =
  let r = ref x in
  for i = 0 to Array.length a - 1 do
    r := f !r i (Array.unsafe_get a i)
  done;
  !r
;;


(*
let fold_lefti f x a = 
  let rec fold_lefti_r xx c = match c with
      0 -> f xx 0 a.(0)
    | _ -> f (fold_lefti_r xx (c-1)) c a.(c)
  in
    fold_lefti_r x ((Array.length a) - 1);;
*)


let array_iter2 f a b =
  let l = Array.length a in
    if l <> Array.length b then invalid_arg "Aux.array_iter2: different lengths"
    else if l <> 0 then 
      begin
	  for i = 0 to l - 1 do
	    f (Array.unsafe_get a i) (Array.unsafe_get b i)
	  done;
      end
;;



let array_map2 f a b =
  let l = Array.length a in
    if l <> Array.length b then invalid_arg "Aux.array_map2: different lengths"
    else if l = 0 then [||] else
      begin
	let r = Array.make l (f (Array.unsafe_get a 0) (Array.unsafe_get b 0) ) in
	  for i = 1 to l - 1 do
	    Array.unsafe_set r i (f (Array.unsafe_get a i) (Array.unsafe_get b i))
	  done;
	  r
      end
;;


let fold_left2 f a b c =
  let l = Array.length a in
    if l <> Array.length b then invalid_arg "Aux.fold_left2: different array lengths"
    else if l = 0 then c else
      begin
	let r = ref c in
	  for i = 0 to l - 1 do
	    r := f (Array.unsafe_get a i) (Array.unsafe_get b i) !r
	  done;
	  !r
      end
;;


let fold_left2i f a b c =
  let l = Array.length a in
    if l <> Array.length b then invalid_arg "Aux.fold_left2: different array lengths"
    else if l = 0 then c else
      begin
	let r = ref c in
	  for i = 0 to l - 1 do
	    r := f i (Array.unsafe_get a i) (Array.unsafe_get b i) !r
	  done;
	  !r
      end
;;




(* splices position x out of array a *)
let splice_pos a y =
  let f z i ai = if ai = y then i else z in
  let x = fold_lefti f (-1) a in
    try
      let b = Array.make ((Array.length a) - 1) 0 in
	Array.blit a 0 b 0 x;
	Array.blit a (x+1) b x ((Array.length a) - x - 1);
	b
    with
	Invalid_argument s -> failwith "splice failed: invalid argument"
      | _ -> failwith "splice failed: problem unknown";;



(* point by point multiplication of two one-dimensional arrays *)
let point_mult a b = array_map2 ( *. ) a b
;;



(* inner product of two one-dimensional arrays *)
let vec_mult a b = 
  let l = Array.length a in
    if l <> Array.length b then invalid_arg "vec_mult : different lengths"
    else if l = 0 then 0. else
      begin
	let r = ref 0. in
	  for i = 0 to l - 1 do
	    r := !r +. a.(i) *. b.(i)
	  done;
	  !r
      end
;;


(* matrix multiplication from the left : y = M x
  a -- M x N array
  b -- N x 1 array
*)
let mat_mult_left a b = 
  let m = Array.length a
  and n = Array.length b in
    if m = 0 || n = 0 then invalid_arg "mat_mult_left : zero dimension "
    else if Array.length a.(0) <> n then invalid_arg "mat_mult_left : dimension don't correspond "
    else
      let c = Array.make m 0. in
	for i = 0 to m-1 do
	  c.(i) <- (vec_mult a.(i) b)
	done;
	c
;;

(* matrix multiplication from the right : y = x M 
   b -- N x M array
   a -- 1 x N array
*)
let mat_mult_right b a = 
  let n = Array.length a in
    if Array.length b <> n then invalid_arg "mat_mult_right : dimension don't correspond "
    else
      let m = Array.length b.(0) in
	if m = 0 || n = 0 then invalid_arg "mat_mult_right : zero dimension "
	else
	  let c = Array.make m 0. in
	  let r = ref 0. in
	    for j = 0 to m-1 do
	      r := 0.;
	      for i = 0 to n - 1 do
		r := !r +. a.(i) *. b.(i).(j)
	      done;
	      c.(j) <- !r;
	    done;
	    c
;;



(* return the list of indices of the elements in array that satisfy a predicate *)
let array_find p arr = 
  let f l i x = if p x then i::l else l in
    List.rev (fold_lefti f [] arr)


(* =============================================== *)



(* =============================================== *)
(* ======== Hash Functions  ==================== *)
(* =============================================== *)

(* map function for hashes *)
let hash_map f h = 
  let h1 = Hashtbl.create (Hashtbl.length h) in
  let func k v =
    Hashtbl.add h1 k (f k v)
  in
    Hashtbl.iter func h;
    h1
;;

let hash_keys h = Hashtbl.fold (fun a b l -> a::l) h []

let hash_vals h = Hashtbl.fold (fun a b l -> b::l) h []

let hash_to_list h = Hashtbl.fold (fun a b l -> (a,b)::l) h []

(* =============================================== *)


(* =============================================== *)
(* ======== Array of Hash and
   ======== array of hash2 functions ============= *)
(* =============================================== *)


(* iter on a hash array*)
let hash_arr_iter f arr_of_hash =
  let get_stat i k v = f i k v in
    Array.iteri (fun i h -> Hashtbl.iter (get_stat i) h) arr_of_hash;;


(* fold on a hash array*)
let hash_arr_fold f arr_of_hash init_val =
  let get_stat i k v s = f i k v s in
    fold_lefti (fun s i h -> Hashtbl.fold (get_stat i) h s) init_val arr_of_hash;;


(* iter on a 2d hash array *)
let hash2_arr_iter f arr_of_2d_hash =
  let get_stat i k1 k2 v = f i k1 k2 v in
    Array.iteri (fun i h -> Hash2.iter (get_stat i) h) arr_of_2d_hash;;


(* fold on a 2d hash array *)
let hash2_arr_fold f arr_of_2d_hash init_val =
  let get_stat i k1 k2 v s = f i k1 k2 v s in
    fold_lefti (fun s i h -> Hash2.fold (get_stat i) h s) init_val arr_of_2d_hash;;


(* =============================================== *)




(* =============================================== *)
(* ======== String Functions  ==================== *)
(* =============================================== *)

(* join string list with a token *)
let string_join token str_l = match str_l with
    [] -> ""
  | [hd] -> hd
  | hd::tl ->
      let m_tl = List.map (function x -> token^x) tl in
	List.fold_left (fun s x -> s^x) hd m_tl;;



(* fold on two sequences*)
let string_fold_i_2 f s1 s2 a =
  if (String.length s1 = String.length s2)
  then
    let rec string_fold_r aa i =
      let c1 = String.sub s1 i 1
      and c2 = String.sub s2 i 1 in
	if (i = String.length s1  - 1)
	then
	  f i c1 c2 aa
	else
	  string_fold_r (f i c1 c2 aa) (i+1)
    in
      string_fold_r a 0
  else
    failwith "string_fold_i_2 failed: Strings have different length!";;

(* =============================================== *)



(* =============================================== *)
(* ======== Mathematics  =========== *)
(* =============================================== *)


let round x = match modf x with
    (x, y) when x < 0.5 -> int_of_float y
  | (x, y) -> (int_of_float y) + 1
;;

let round_f x = match modf x with
    (x, y) when x < 0.5 -> y
  | (x, y) ->  y +. 1.
;;


let loga a x = match a, x with
    0., _ -> nan
  | 1., _ -> nan
  | _, 0. -> neg_infinity
  | _, x when x < 0. -> nan
  | _, _ -> (log x) /. (log a)
;;


(* product of natural numbers between n1 and n2 inclusive *)
let prod_nat n1 n2 =
  let (m1, m2) = match (n1, n2) with
      (x, y) when x < 0 || y < 0 -> failwith "prod_nat: negative argument(s)"
    | (0, 0) -> (1, 1)
    | (x, y) when x > y -> (1, 1)
    | (0, y) -> (1, y)
    | _ -> (n1, n2)
  in
  let rec multiply p n =
    if n = m2
    then Big_int.mult_int_big_int n p
    else multiply (Big_int.mult_int_big_int n p) (n + 1)
  in
    multiply (Big_int.big_int_of_int 1) m1;;


(* factorial of natural numbers *)
let fact = prod_nat 1;;


(* n choose x *)
let choose n x =
  if n >= x then
    let y = n - x in
      if y >= x then
	Big_int.float_of_big_int (Big_int.div_big_int (prod_nat (y+1) n) (fact x))
      else
	Big_int.float_of_big_int (Big_int.div_big_int (prod_nat (x+1) n) (fact y))
  else
    failwith "get_choose: n < x";;



(* binomial probability *)
let binom_p n p x = (choose n x) *. p**(float_of_int x) *. (1. -. p)**(float_of_int (n - x));;


(* binomail cumulative distribution function *)
let binom_cdf n p x =
  let l = Array.to_list ( Array.init ((int_of_float x) + 1) ((+) 0) ) in
    List.fold_left (fun s y -> s +. binom_p n p y) 0. l;;


(* generates a list of m different random integers between 0 (incl) and n (excl) *)
let gen_rand_l n m =
  if m <= n then
    let rec gen arr l rem = match rem with
	0 -> l
      | _ ->
	  let b = Array.length arr in
	  let i = Random.int b in
	  let new_arr = Array.append (Array.sub arr 0 i) (Array.sub arr (i+1) (b-i-1)) in 
	    gen new_arr (i::l) (rem-1)
    in
      gen (Array.init n ((+) 0) ) [] m
  else
    failwith "gen_rand_l : m > n ! \n"
;;
    


(* generate a number from the Poisson distribution (Knuth) *)
let gen_poisson lambda =
  if lambda > 0. then
    let r = exp (-. lambda)
    and p = ref 1.
    and k = ref 0 in
      while !p > r do
	k := !k + 1;
	p := !p *. (Random.float 1.)
      done;
      !k - 1
  else
    failwith "gen_poisson: negative rate\n"
;;



(* generate a number from the exponential distribution *)
let gen_exp a =
  let x = Random.float 1. in
    -. (log x) /. a
;;



(* === generate a number from a discrete distribution === *)

let print_distrib l =
  let f (cp, i) = print_string ((string_of_float cp)^" "^(string_of_int i)^"\n") in
    List.iter f l
;;


let get_distrib arr =
  let f (r, l) i p =
    let new_r = r +. p in
      (new_r, (new_r, i)::l)
  in
  let (final_r, final_l) = fold_lefti f (0., []) arr in
    (final_r, List.rev final_l)
;;



let rec draw distrib x =
  (* print_string ((string_of_float x)^" "); *)
  match distrib with
      [] -> failwith ("draw: failing x "^(string_of_float x)^"the random number fell out of range")
    | (cp, i)::[] -> i
    | (cp, i)::tl -> if x < cp then i else (draw tl x)
;;



let gen_discr arr =
  let (total_r, distrib) = get_distrib arr  in
    draw distrib (Random.float total_r)
;;






let gen_rand_float_l r n =
  let rec gen i l = match i with
      0 -> l
    | _ -> gen (i-1) ((Random.float r)::l)
  in
    gen n []
;;


let draw_n distrib rand_list =
  let rec draw_n_r distr rand_l outp_l =
    match rand_l with
	[] -> outp_l
      | r_hd::r_tl -> match distr with
	    [] -> failwith ("draw_n: failing to find number "^(string_of_float r_hd)^" in the range")
	  | (cp, i)::tl -> if r_hd < cp then (draw_n_r distr r_tl (i::outp_l)) else (draw_n_r tl rand_l outp_l)
  in
    draw_n_r distrib rand_list []
;;




let gen_discr_n arr n =
  let (total_r, distrib) = get_distrib arr  in
  let rand_list = List.sort compare (gen_rand_float_l total_r n) in
(*    print_endline (String.concat " " (List.map string_of_float rand_list));
    print_float total_r; *)
    List.rev (draw_n distrib rand_list)
;;


(* mean value of a list *)
let mean = function 
    [] -> raise (Invalid_argument "mean : empty list")
  | l -> let s = List.fold_left ( +. ) 0. l in s /. (float_of_int (List.length l))
;;
    


(* =============================================== *)
