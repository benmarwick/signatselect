(* checks if a keys k1 and k2 are bound in the 2D hash h *)
let mem h k1 k2 =
  if (Hashtbl.mem h k1)
  then
    Hashtbl.mem (Hashtbl.find h k1) k2
  else
    let () = Hashtbl.add h k1 (Hashtbl.create 1) in
      false;;

(* get a value from a 2D hash *)
let find h k1 k2 =
  Hashtbl.find (Hashtbl.find h k1) k2;;


(* set a value in a 2D hash *)
let replace h k1 k2 v =
  if not (Hashtbl.mem h k1) then Hashtbl.add h k1 (Hashtbl.create 1);
  Hashtbl.replace (Hashtbl.find h k1) k2 v;;

(* iter on a 2D hash *)
let iter f h =
  let myf a b =
    Hashtbl.iter (f a) b
  in
    Hashtbl.iter myf h;;

(* folds a 2D hash *)
let fold f h x =
  let myf a b y =
    Hashtbl.fold (f a) b y
  in
    Hashtbl.fold myf h x;;





let keys h =
  let f a b c l = (a,b)::l in
    fold f h []
  


