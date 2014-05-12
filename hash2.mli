(** This module contains definitions and functions for two-dimensional hash tables *)

(** {3 General functions} *)

val keys : ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> ('a * 'b) list 
(** [keys h] returns the list of keys of the hash [h] *)

(** {3 Iterators }*)

val mem : ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> 'a -> 'b -> bool
(** [hash2_mem h k1 k2] checks if the value of the two-dimensional hash [h] corresponding to keys [k1] and [k2] exists. If key [k1] is not bound, binds it with a new hash *)

val find : ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> 'a -> 'b -> 'c;;
(** [hash2_find h k1 k2] returns the value of the two-dimensional hash [h] corresponding to keys [k1] and [k2]. It is equivalent to [ Hashtbl.find (Hashtbl.find h k1) k2] *)

val replace : ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> 'a -> 'b -> 'c -> unit
(** [hash2_set h k1 k2 v] replaces the existing value of the two-dimensional hash [h] corresponding to keys [k1] and [k2] with [v] or, if this entry does not yet exist, creates it. If the key [k1] is already bound, this function is equivalent to [Hashtbl.replace (Hashtbl.find h k1) k2 v] *) 

val iter : ('a -> 'b -> 'c -> unit) -> ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> unit;;
(** [hash_iter_2 f h] applies function [f k1 k2 v] to each element of hash [h]. Here, [v] corresponds to the binding of [k1 k2] *)

val fold : ('a -> 'b -> 'c -> 'd -> 'd) -> ('a, ('b, 'c) Hashtbl.t) Hashtbl.t -> 'd -> 'd;;
(** [hash_fold_2 f h x] folds a 2d hash with function [f k1 k2 v y] where [v] corresponds to the binding of [k1 k2]. [x] is the initial fold value *)
