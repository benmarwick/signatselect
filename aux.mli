(** This module encompasses various useful auxilary functions *)

(** {2 List Functions} *)

val splice : 'a list -> int -> int -> 'a list
(** [splice l o len] returns the list of length [len] consisting of consecuitve elements of list [l] starting from the element at position [o] of [l] *)

val purge : 'a list -> 'a list
(** [purge l] removes all duplicate entries from a sorted list [l] (i.e., all duplicate entries must be adjacent) *)


(** {2 Array Functions}*)

val matrix_tr : 'a array array -> 'a array array;;
(** [matrix_tr m ] returnes the transpose of matrix [m]. *)

val matrix_copy : 'a array array -> 'a array array;;
(** [matrix_copy m ] returnes a fresh copy of the matrix [m]. The elements of the resulting matrix are completely unaffected by changes in [m], and vise versa. Note that this would not be the case if one simply did [Array.copy m].*)

val fold_lefti : ('a -> int -> 'b -> 'a) -> 'a -> 'b array -> 'a;;
(** [fold_lefti f x a] folds an array with index. It is equivalent to [ f ( ... (f (f x 0 a.(0)) 1 a.(1) ) ...) n a.(n) ] *)

val array_iter2 : ('a -> 'b -> unit) -> 'a array -> 'b array -> unit ;;
(** [array_iter2 f a1 a2] iterates through two arrays. It is equivalent to [f a1.(0) a2.(0);  f a1.(1) a2.(1); ... ; f a1.(n) a2.(n)]. Raises [Invalid_argument] if [a1] and [a2] have different lengths. *)


val array_map2 : ('a -> 'b -> 'c) -> 'a array -> 'b array -> 'c array;;
(** [array_map2 f a1 a2] maps two arrays with a function. It is equivalent to [[| f a1.(0) a2.(0);  f a1.(1) a2.(1); ... ; f a1.(n) a2.(n) |]]. Raises [Invalid_argument] if [a1] and [a2] have different lengths. *)


val fold_left2 : ('a -> 'b -> 'c -> 'c) -> 'a array -> 'b array -> 'c -> 'c;;
(** [fold_left2 f a1 a2 b] folds two arrays with a function. It is equivalent to [ f a1.(n) a2.(n) ( ... (f a1.(0) a2.(0) b) ...)]. Raises [Invalid_argument] if [a1] and [a2] have different lengths. *)

val fold_left2i : (int -> 'a -> 'b -> 'c -> 'c) -> 'a array -> 'b array -> 'c -> 'c;;
(** [fold_left2 f a1 a2 b] folds two arrays with index with a function. It is equivalent to [ f n a1.(n) a2.(n) ( ... (f 0 a1.(0) a2.(0) b) ...)]. Raises [Invalid_argument] if [a1] and [a2] have different lengths. *)

val array_find : ('a -> bool) -> 'a array -> int list ;;
(** [array_find p arr] returns the list of indices of elements of array [arr] that satisfy the predicate [p]. Indices are ordered from lowest to highest*)


val splice_pos : int array -> int -> int array;;
(** [splice_pos a y] splices out of the array the last position in the array that contains element y *)

val point_mult : float array -> float array -> float array ;;
(** [point_mult a b] returns the element-by-element product of two one-dimensional arrays. Arrays must have the same length *)

val vec_mult : float array -> float array -> float ;;
(** [vec_mult a b] returns the inner product of two one-dimensional arrays. Arrays must have the same length *)

val mat_mult_left : float array array -> float array -> float array ;;
(** [mat_mult m v] multiplies matrix [m] by the vector [v] from the left : [m x b]. [m] is an [m x n] matrix, [v] is an array of length [n]. *)

val mat_mult_right : float array array -> float array -> float array ;;
(** [mat_mult m v] multiplies matrix [m] by the vector [v] from the right : [v m]. [m] is an [n x m] matrix, [v] is an array of length [n]. *)




(** {2 Hash Functions}*)

val hash_map : ('a -> 'b -> 'c) -> ('a, 'b) Hashtbl.t -> ('a, 'c) Hashtbl.t
(** [hash_map f h] maps hash [h] into a new hash with the same keys and values determined by [f k v] where [k] and [v] are the keys and values of [h], respectively *)

val hash_keys : ('a, 'b) Hashtbl.t -> 'a list
(** [hash_keys h] returns the list of hash keys *)

val hash_vals : ('a, 'b) Hashtbl.t -> 'b list
(** [hash_keys h] returns the list of hash values *)

val hash_to_list : ('a, 'b) Hashtbl.t -> ('a * 'b) list
(** [hash_keys h] returns the list pairs of keys and values *)



(** {2 Hash Array and Hash2 Array Functions}*)

val hash_arr_iter : (int -> 'a -> 'b -> unit) -> ('a, 'b) Hashtbl.t array -> unit
(** [hash_arr_iter f ha] iterates on a hash array; [f i k v] where [i] -- array index, [k] -- hash key, [v] -- hash value *)

val hash_arr_fold : (int -> 'a -> 'b -> 'c -> 'c) -> ('a, 'b) Hashtbl.t array -> 'c -> 'c
(** [hash_arr_fold f ha c] folds a hash array; [f i k v s] where [i] -- array index, [k] -- hash key, [v] -- hash value, [s] -- fold value; [c] -- initial fold value *)

val hash2_arr_iter : (int -> 'a -> 'b -> 'c -> unit) -> ('a, ('b, 'c) Hashtbl.t) Hashtbl.t array -> unit
(** [hash2_arr_iter f h2a] iterates on a hash2 array; [f i k1 k2 v] where [i] -- array index, [k1] and [k2] -- hash keys, [v] -- hash value *)

val hash2_arr_fold : (int -> 'a -> 'b -> 'c -> 'd -> 'd) -> ('a, ('b, 'c) Hashtbl.t) Hashtbl.t array -> 'd -> 'd
(** [hash2_arr_fold f h2a c] folds a hash2 array; [f i k1 k2 v s] where [i] -- array index, [k1] and [k2] -- hash keys, [v] -- hash value, [s] -- current fold value, [c] -- initial fold value *)






(** {2 String Functions}*)

val string_join : string -> string list -> string
(** [string_join t [s1; s2; ...; sn]] join the string list with a token [t]. It is equivalent to [s1^s2^...^sn] *)

val string_fold_i_2 : (int -> string -> string -> 'a -> 'a) -> string -> string -> 'a -> 'a;;
(** [string_fold_i_2 f s1 s2 x] folds two strings of equal length with index. It is equivalent to [ f n (String.sub s1 n 1) (String.sub s2 n 1) ( ... (f 1 (String.sub s1 1 1) (String.sub s2 1 1) (f 0 (String.sub s1 0 1) (String.sub s2 0 1) x)) ...)  ]. Note that individual characters of the string are extracted as strings, not chars. *)








(** {2 Mathematics } *)

val round : float -> int
(** [round x] rounds [x] to the nearest integer *)

val round_f : float -> float
(** [round x] rounds [x] to the nearest integer but keeps type [float]*)

val loga : float -> float -> float
(** [loga a x] : logarithm x base a *)

val prod_nat : int -> int -> Big_int.big_int
(** [prod_nat n1 n2] computes the product of integers between [n1] and [n2] inclusive. If one or both arguments are negative, the function fails. The following convetions are used:
 1) [prod_nat 0 0 = 1]
 2) [prod_nat n1 n2 = 1] if [n1 > n2] 
 3) [prod_nat 0 n2 = prod_nat 1 n2] *)

val fact : int -> Big_int.big_int
(** [fact n] computes the factorial for n. By convention [fact 0 = 1] *)

val choose : int -> int -> float
(** [choose n x] computes "n choose x", i.e., n!/x!/(n-x)! *)

val binom_p : int -> float -> int -> float
(** [binom_p n p x] computes the probability of [x] successes in a series of [n] binomial trials with success probability [p] *)

val binom_cdf : int -> float -> float -> float 
(** [binom_cdf n p x] computes the cumulative distribution function for the binomial distribution with parameters [n] and [p], at point [x] *)

val gen_rand_l : int -> int -> int list
(** [gen_rand_l n m] generates a list of [m] different random numbers between [0] and [n].[m] must not exceed [n] *)

val gen_poisson : float -> int
(** [gen_poisson lambda] generates a random number from the Poisson distribution with parameter [lambda]. Uses Knuth's algorithm. *)

val gen_exp : float -> float
(** [gen_exp a] generates a random number from an exponential distribution with parameter [a] *)

val gen_discr : float array -> int
(** [gen_discr p_vec] generates a number from a discrete distribution specified by the probabaility distribution vector [p_vec]. Note that [p_vec] DO NOT need to sum up to 1. Value [i] is drawn with weight [p_vec.(i)] *)

val gen_discr_n : float array -> int -> int list
(** [gen_discr_n p_vec n] generates a list of [n] numbers from a discrete distribution specified by the probabaility distribution vector [p_vec].  Note that [p_vec] DO NOT need to sum up to 1. Value [i] is drawn with weight [p_vec.(i)] *)

val mean : float list -> float
(** [mean l] computes the mean value of the list. Raises [Invalid_argument "mean : empty list"] if list is empty *)
