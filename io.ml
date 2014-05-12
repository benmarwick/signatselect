open Printf;;

exception Invalid_format of string;; (* Exception raised if the file format is invalid *)


(*===================================================*)
(*============= INPUT operations =====================*)
(*===================================================*)

(* reading the file _filename_ into a list of lines
   in: filename
   out: list of lines (type string) withou new line character *)
let read_file_to_list filename = 
  try 
    let infile = open_in filename
    and a = ref []
    in
      begin
	try
	  while true do
	    a := !a@[(input_line infile)];
	  done
	with End_of_file -> ();
      end;
      close_in infile;
      !a;
  with
      _ -> let () = print_string "Could not read the file" in [];;

(* parses the fasta file:*)
let get_seqs_fasta filename =
  let input = open_in filename in
  let s = Stream.of_channel input in
  let h = Hashtbl.create 1
  and id = ref ""
  and seq = ref ""
  and read_id = ref false
  and read_seq = ref false
  and eol = ref false in
  let reset_all () =
    if !id <> ""
    then
      let my_id::tl = Str.split (Str.regexp "/") !id in
	Hashtbl.add h my_id !seq;
	id := ""; seq := ""; read_id := false; read_seq := false; eol := false
  in
  let myf c = match (c, !read_id, !read_seq) with
      (* --- reading nothing --- *)
      ('>', false, false) -> read_id := true; eol := false;
    | (_, false, false) -> ()

      (* --- reading sequence --- *)
    | ('\n', false, true) -> if !eol then reset_all () else eol := true
    | ('>', false, true) -> reset_all (); read_id := true
    | (_, false, true) -> eol := false; seq := !seq^(Char.escaped (Char.uppercase c))

      (* --- reading id --- *)
    | ('\n', true, false) ->
	read_id := false; read_seq := true
    | (_, true, false) -> id := !id^(Char.escaped c)

      (* --- errors --- *)
    | (_, true, true) -> failwith "cannot read sequence and id simultaneously!"
  in
    Stream.iter myf s;
    reset_all ();
    close_in input;
    h;;


(*let rec stream_fold f s x = match s with parser
    [< 'y >] -> stream_fold f s (f y x)
  | [< >] -> x;;
*)


(* parses the fasta file:*)
(* let get_seqs_fasta filename =
  let name_match = Str.regexp "^>\\([^/]+\\)"
  and file = read_file_to_list filename in
  let h = Hashtbl.create (List.length file) in
  let rec f id seq line_list = match line_list with
      [] -> ()
    | line::tl ->
	(* print_string (id^" "^seq^"\n"); *)
	match id with
	  "" ->
	    if Str.string_match name_match line 0
	    then f (Str.matched_group 1 line) seq tl
	    else f "" "" tl
	| _ ->
	    if (line = "") || (Str.string_match (Str.regexp "^ +") line 0)
	    then
	      match seq with
		  "" -> failwith ("Empty sequence "^id^"\n")
		| _ ->
		    Hashtbl.add h id seq;
		    f "" "" tl
	    else
(*	      if Str.string_match (Str.regexp "[^ARNDCEQGHILKMFPSTWYVarndceqghilkmfpstwyv]") line 0
	      then *)
	      if (Str.string_match name_match line 0)
	    then
	      begin
		Hashtbl.add h id seq;
		f (Str.matched_group 1 line) "" tl
		end
(*	    else
		begin
		  print_string (line^"\n");
		  failwith "Invalid FASTA format: encountered an unknown character!"
	    end *)
	    else
	      f id (seq^(String.uppercase line)) tl
  in
    f "" "" file;
    h;;
*)




(* parses the interaction file into the actual interaction matrix:
   in: filename
   out: hash of hashes of interactions *)
let get_j filename =
  try  
    let list_of_lines = read_file_to_list filename in
    let jj = Hashtbl.create (List.length list_of_lines) in
    let assign_j line =
      let name::neighb = Str.split (Str.regexp ":") line in
	match neighb with
	    [] -> Hashtbl.add jj name (Hashtbl.create 0)
	  | [neighbors] ->
	      let neighb_j = Str.split (Str.regexp " ") neighbors in
	      let rec iterate l key  = match l with
		  [] -> ()
		| hd::tl when (key = "") -> iterate tl hd
		| hd::tl -> Hashtbl.add (Hashtbl.find jj name) key (float_of_string hd); iterate tl ""
	      in
		Hashtbl.add jj name (Hashtbl.create 1);
		iterate neighb_j "";
    in
      List.iter assign_j list_of_lines;
      jj
  with
      _ ->
	print_string "Io.get_j: Could not read the file!\n";
	exit 1
;;



(* parses the metric file into the actual metric:
   in: filename
   out: hash of hashes of distances *)
let get_metric filename =
  let title::triag = read_file_to_list filename in
  let names = Str.split (Str.regexp "\t") title in
  let mat = Hashtbl.create (List.length names) in
  let assign_matrix distance_string name1 =
    let rec iterate dist_list name_list = match (dist_list, name_list) with
	([], l) -> ();
      | (d::tl1, name2::tl2) ->
	  Hashtbl.add (Hashtbl.find mat name1) name2 (float_of_string d);
	  Hashtbl.add (Hashtbl.find mat name2) name1 (float_of_string d);
	  iterate tl1 tl2;
    and distances_x = Str.split (Str.regexp "\t") distance_string
    in
      Hashtbl.add mat name1 (Hashtbl.create (List.length names));
      iterate distances_x names;
  in
    List.iter2 assign_matrix triag names;
    mat;;







(* parses the input line into parameters and options *)
let parse_input_line line_arr =
  let line_list = Array.to_list line_arr
  and get_params_options l =
    let rec get_params_r all params = match all with
	[] -> (params, all)
      | x::tl -> match (String.get x 0) with
	    '-' ->
	      let options = List.map (function s -> String.sub s 1 (String.length s - 1)) all in
		(params, options)
	  | _ -> get_params_r tl (params@[x])
    in
      match l with
	  [] -> failwith "empty input line!";
	| hd::tl -> get_params_r tl []
  in
    get_params_options line_list;;
    
      

(* parses a parameter file: each parameter name must be marked with an '>'; to each parameter then corresponds an list of lines until the first line gap *)
let parse_param_file filename =
  if Sys.file_exists filename
  then
    let line_list = read_file_to_list filename in
    let params = Hashtbl.create 1
    and finished_params = Hashtbl.create 1 in
    let check_line id line = match (id, line) with
	("", "") -> ""
      | ("", _) ->
	  if String.get line 0  = '>'
	  then
	    let new_id = String.sub line 1 ((String.length line) - 1) in
	      if new_id = "" then failwith "empty id!";
	      new_id
	  else
	    ""
      | (_, _) ->
	  if Str.string_match (Str.regexp "^ *$") line 0
	  then
	    begin
	      Hashtbl.add finished_params id 1;
	      ""
	    end
	  else
	    begin
	      if String.get line 0  = '>'
	      then
		let new_id = String.sub line 1 ((String.length line) - 1) in
		  Hashtbl.add finished_params id 1;
		  if new_id = "" then failwith "empty id!";
		  new_id
	      else
		begin
		  if Hashtbl.mem params id
		  then
		    let curr_entry = Hashtbl.find params id in
		      if Hashtbl.mem finished_params id then failwith ("Duplicate parameter entry: "^id);
		      Hashtbl.replace params id (curr_entry@[line])
		  else
		    Hashtbl.add params id [line];
		  id
		end
	    end
    in
    let x = List.fold_left check_line "" line_list in
      params
  else
    failwith ("File "^filename^" does not exist");;




(*
let get_tree_sim_params filename =
  let l = read_file_to_list filename
  and r_m = Hashtbl.create 4 in
  let f (omega, init_seq, output, stop, flags) a = match flags with
      "", false, false, false, false ->
	if  Str.string_match (Str.regexp "^rate_matrix") a 0
	then (omega, init_seq, output, stop, ("A", false, false, false, false))
	else if Str.string_match (Str.regexp "^dnds") a 0
	then (omega, init_seq, output, stop, ("", true, false, false, false))
	else if Str.string_match (Str.regexp "^init_seq") a 0
	then (omega, init_seq, output, stop, ("", false, true, false, false))
	else if Str.string_match (Str.regexp "^output") a 0
	then (omega, init_seq, output, stop, ("", false, false, true, false))
	else if Str.string_match (Str.regexp "^stop") a 0
	then (omega, init_seq, output, stop, ("", false, false, false, true))
	else (omega, init_seq, output, stop, ("", false, false, false, false))
    | c, false, false, false, false when String.length c > 0 ->
	let (l1, next) = match c with
	    "A" -> (['C'; 'G'; 'T'], "C")
	  | "C" -> (['A'; 'G'; 'T'], "G")
	  | "G" -> (['A'; 'C'; 'T'], "T")
	  | "T" -> (['A'; 'C'; 'G'], "")
	and l2 = List.map (function s -> float_of_string s) (Str.split (Str.regexp "\t") a) in
	  Hashtbl.add r_m (String.get c 0) (List.combine l1 l2);
	  (omega, init_seq, output, stop, (next, false, false, false, false))
    | "", true, false, false, false ->
	let l_string = Str.split (Str.regexp "\t") a in
	let a = Array.of_list (List.map (function s -> float_of_string s) l_string) in
	  (a, init_seq, output, stop, ("", false, false, false, false))
    | "", false, true, false, false ->
	(omega, (String.uppercase a), output, stop, ("", false, false, false, false))
    | "", false, false, true, false ->
	(omega, init_seq, a, stop, ("", false, false, false, false))
    | "", false, false, false, true -> (omega, init_seq, output, (int_of_string a), ("", false, false, false, false))
    | _ -> failwith "Wrong file format" in
  let (my_dnds, my_init_seq, my_output, my_stop, my_flags) = List.fold_left f ([||], "", "", 1, ("", false, false, false, false)) l in
    (r_m, my_dnds, my_init_seq, my_output, my_stop);;
*)

(* ============= done =========================== *)


(* ================================================ *)
(* ============== OUTPUT FUNCTIONS ==================*)
(* ================================================ *)

let print_metric_f keys h filename =
  let out = open_out filename in
  let output_l l =
    output_string out (sprintf "%s\n" (String.concat "\t" l));    
  in
    output_l keys;
    List.iter ( function k -> output_l (List.map (sprintf "%.5g") (Hashtbl.find h k))  ) keys;
    close_out out
;;





let print_metrics (keys, h) filename =
  let cod_out = open_out (filename^".cod.met")
  and aa_out = open_out (filename^".aa.met")
  and sil_out = open_out (filename^".sil.met")
  and make_printable s = match s with
      "" -> "\n"
    | s -> (String.sub s 1 ( (String.length s)-1 ))^"\n"
  in
  let line1 = make_printable (List.fold_left (fun s x -> s^"\t"^x ) ""  keys)
  and print_numbers key =
    let l = Hashtbl.find h key in
    let f (cod_s, aa_s, sil_s) (cod_d, aa_d, sil_d) =
      (cod_s^"\t"^(string_of_int cod_d), aa_s^"\t"^(string_of_int aa_d), sil_s^"\t"^(string_of_int sil_d)) in
    let (cod_string, aa_string, sil_string) = List.fold_left f ("", "", "") l in
      output_string cod_out (make_printable cod_string);
      output_string aa_out (make_printable aa_string);
      output_string sil_out (make_printable sil_string)
  in
    output_string cod_out line1;
    output_string aa_out line1;
    output_string sil_out line1;
    List.iter print_numbers keys;
    close_out cod_out;
    close_out aa_out;
    close_out sil_out;;
    

let print_metrics_f (keys, h) filename =
  let cod_out = open_out (filename^"cod.met")
  and aa_out = open_out (filename^"aa.met")
  and sil_out = open_out (filename^"sil.met")
  and make_printable s = match s with
      "" -> "\n"
    | s -> (String.sub s 1 ( (String.length s)-1 ))^"\n"
  in
  let line1 = make_printable (List.fold_left (fun s x -> s^"\t"^x ) ""  keys)
  and print_numbers key =
    let l = Hashtbl.find h key in
    let f (cod_s, aa_s, sil_s) (cod_d, aa_d, sil_d) =
      (cod_s^"\t"^(string_of_float cod_d), aa_s^"\t"^(string_of_float aa_d), sil_s^"\t"^(string_of_float sil_d)) in
    let (cod_string, aa_string, sil_string) = List.fold_left f ("", "", "") l in
      output_string cod_out (make_printable cod_string);
      output_string aa_out (make_printable aa_string);
      output_string sil_out (make_printable sil_string)
  in
    output_string cod_out line1;
    output_string aa_out line1;
    output_string sil_out line1;
    List.iter print_numbers keys;
    close_out cod_out;
    close_out aa_out;
    close_out sil_out;;



let print_file_hash2 h filename =
  let outfile = open_out filename in
    Hashtbl.iter (
      fun a b ->
	output_string outfile (a^":");
	Hashtbl.iter (fun c d -> output_string outfile (c^" "^(string_of_float d)^" ") ) b;
	output_string outfile "\n";
    ) h;
    close_out outfile;;



let print_int_hash_file h filename =
  let output = open_out filename in
    Hashtbl.iter (fun a b -> output_string output (a^" "^(string_of_int b)^"\n") )  h;
    close_out output
;;



let print_int_hash h =
  Hashtbl.iter (fun a b -> print_string (a^" "^(string_of_int b)^"\n") )  h;;

let print_float_hash h = Hashtbl.iter (fun a b -> print_string (a^" "^(string_of_float b)^"\n") ) h;;

let print_string_hash h = Hashtbl.iter (fun a b -> print_string (a^" "^b^"\n")) h;;


let print_float_hash2 h =
  Hashtbl.iter (
    fun a b ->
      print_string (a^" ");
      Hashtbl.iter (fun c d -> print_string (c^" "^(string_of_float d)^" ") ) b;
      print_string "\n";
  ) h;;


let print_seqs_fasta h filename =
  let outfile = open_out filename in
  let f a b = output_string outfile (">"^a^"\n"^b^"\n\n") in
    Hashtbl.iter f h;
    close_out outfile;;





(* ================ done ============================*)
