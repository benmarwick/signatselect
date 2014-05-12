(** This module implements Input/Output operations useful for genetic sequence analysis, Super-Paramagnetic Clustering and phylogenetic analysis *)

(** {2 Exceptions }*)

exception Invalid_format of string;;
(** Exception raised if the file format is invalid *)


(** {2 Input methods} *)

val read_file_to_list : string -> string list;;
(** [read_file_to_list filename] reads file [filename] to a list of lines chopping off the newline characters *)

val get_seqs_fasta : string -> (string, string) Hashtbl.t;;
(** [get_seqs_fasta filename] reads the file [filename] in the FASTA format for both aminoacid and nucleotide sequences and returns  a hash of sequences *)

val get_j : string -> (string, (string, float) Hashtbl.t) Hashtbl.t;;
(** [get_j filename] reads file [filename] containing the neighborhoods of points into the hash of neighborhoods  *)

val get_metric : string -> (string, (string, float) Hashtbl.t) Hashtbl.t;;
(** [get_metric filename] reads file [filename] containing the metric on some objects into the hash of hashes of distances  *)


val parse_input_line : string array -> string list * string list;;
(** [parse_input_line input_line_array] parses the input line (obtained by [Sys.argv]) into input parameters and input options. Options are defined as values of [input_line_array] that start with a '-'). Input parameters are defined as all values of the [input_line_array] between 1 (the name of the program is dropped) and the first input option. For example, the line [my_program 12 8 -N:100 -f:extra] would return [(["12"; "8"], ["N:100"; "f:extra"]) ]. *)

val parse_param_file : string -> (string, string list) Hashtbl.t;;
(** [parse_param_file line_list] parses the line_list so that each line marked with a '>' becomes a key (=parameter name). All lines after the parameter name before the first gap line become the value in the resulting hash table in the entry corresponding to the paramter name *)

(*
val get_tree_sim_params : string -> (char, (char * float) list) Hashtbl.t * float array * string * string * int;;
(** [get_tree_sim_params filename] reads the parameters for the simulation of the evolution along a phylogenetic tree by the program e_sim. The output is a 5-tupole: [(rate_matrix, dnds_array, initial_sequence, output_filename, stop)]. Here [stop] is a variable that indicated whether or not the stop codons should be treated as such ([stop <> 0]) or not ([stop = 0]) *)
*)



(** {2 Output methods} *)

val print_int_hash : (string, int) Hashtbl.t -> unit;;
(** [print_int_hash hash] prints the [hash] to the standard output *)

val print_float_hash : (string, float) Hashtbl.t -> unit;;
(** [print_int_hash hash] prints the [hash] to the standard output *)

val print_string_hash : (string, string) Hashtbl.t -> unit;;
(** [print_int_hash hash] prints the [hash] to the standard output *)

val print_float_hash2 : (string, (string, float) Hashtbl.t) Hashtbl.t -> unit;;
(** [print_float_hash2 hash] prints the 2-dimensional [hash] to the standard output*)

val print_int_hash_file : (string, int) Hashtbl.t -> string -> unit;;
(** [print_int_hash hash] prints the [hash] to the file [filename] *)

val print_metric_f : string list -> (string, float list) Hashtbl.t -> string -> unit;;
(** [print_metric_f keys metric filename] prints the [metric] into the file [filename] *)


val print_metrics : string list * (string, (int*int*int) list) Hashtbl.t -> string -> unit;;
(** [print_metrics (keys, metrics) filename] prints the three metrics into the three corresponding files: [filename.cod.met], [filename.aa.met], [filename.sil.met]. [keys] is the list of keys of hash [metrics]. Each binding of [metrics] is a list of 3-tupoles [(codon_dist, aa_dist, sil_dist)] of distances to the members of the hash that stand before the current key in the [keys] list. *)

val print_metrics_f : string list * (string, (float*float*float) list) Hashtbl.t -> string -> unit;;
(** Same as [print_metrics] but with float metric values*)
 
val print_file_hash2 : (string, (string, float) Hashtbl.t) Hashtbl.t -> string -> unit;;
(** [print_file_hash2 hash filename] opens a file with the name [filename] and prints the 2-dimensional [hash] into it in the neighborhood format *)

val print_seqs_fasta : (string, string) Hashtbl.t -> string -> unit;;
(** [print_seqs_fasta h filename] stores a hash of sequences in a file in the FASTA format *)
