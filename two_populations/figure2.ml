open Printf
open Owl

(* code to reproduce Figure 2 *)

module D = Dynamics.Make (Prms.Default)


(* --------------------------------------------------------------------------------
   ---    panel C : transient with input steps                                  ---
   -------------------------------------------------------------------------------- *)

let input = 
  let f t = Mat.create 2 1 (if t < 2. then 0. else if t < 6. then 2. else 15.) in
  D.(Time_varying f)

let us, rs = D.simulate ~verbose:true ~duration:10. ~input ~u0:(Mat.zeros 2 1) ()
let us = Arr.squeeze us and rs = Arr.squeeze rs
let _ = Mat.save_txt (Mat.transpose us) (Dir.in_dir "us_input_steps")
let _ = Mat.save_txt (Mat.transpose rs) (Dir.in_dir "rs_input_steps")


(* --------------------------------------------------------------------------------
   ---    panel D : first- and second-order moments                             ---
   -------------------------------------------------------------------------------- *)

let n_trials = 5000
let burn_in = 1.0
let time_slice = [- Maths.(round (0.1 /. D.sampling_dt) |> int_of_float); -1]
let prepare v = 
  v 
  |> Arr.get_slice [[]; []; time_slice]
  |> Arr.transpose ~axis:[| 1; 0; 2 |] (* neurons first *)
  |> (fun v -> (* collapse last 2 dims *)
      let d = Arr.shape v in Arr.reshape v [| d.(0); d.(1) * d.(2) |])
  |> Mat.transpose

let inputs = Array.init 20 (fun i -> D.(Steady (Mat.create 2 1 (float i))))
let u0 = Mat.zeros 2 n_trials

type result = { u_mean: Mat.mat; u_var: Mat.mat; 
                r_mean: Mat.mat; r_var: Mat.mat }

let results = Array.mapi (fun i input ->
    Gc.compact ();
    printf "input: %i\n%!" i;
    let us, rs = D.simulate ~verbose:true ~duration:burn_in ~input ~u0 () in
    let us = prepare us and rs = prepare rs in
    { u_mean = Mat.mean ~axis:0 us;
      r_mean = Mat.mean ~axis:0 rs;
      u_var = Mat.var ~axis:0 us;
      r_var = Mat.var ~axis:0 rs }
  ) inputs

(* save the results *) 
let _ =
  let save extract filename = 
    results 
    |> Array.map extract 
    |> Mat.concatenate ~axis:0
    |> (fun m -> Mat.(save_txt m (Dir.in_dir filename)))
  in
  save (fun x -> x.u_mean) "u_mean";
  save (fun x -> x.r_mean) "r_mean";
  save (fun x -> x.u_var) "u_var";
  save (fun x -> x.r_var) "r_var"


