open Printf
open Owl

(* code to reproduce Figure 4 *)

module P = Prms.Default
module D = Dynamics.Make (P)

let stim =
  let baseline = 2. in
  let a_max = 20. in
  fun contrast ->
    let z = contrast *. a_max in
    Mat.(baseline $+ (z $* P.stimulus)) 


let () = 
  Mat.save_txt D.w (Dir.in_dir "w");
  Mat.save_txt D.sigma_noise (Dir.in_dir "sigma_noise");
  Mat.save_txt P.stimulus (Dir.in_dir "stimulus")

(* --------------------------------------------------------------------------------
   ---    panel E : transient at stimulus onset                                 ---
   -------------------------------------------------------------------------------- *)

let transient_input t_switch = 
  let low_stim = stim 0. and high_stim = stim 1. in
  let f t = if t<t_switch then low_stim else high_stim in
  D.(Time_varying f)

let us, rs = D.simulate ~verbose:true ~duration:10. ~input:(transient_input 5.) ~u0:(Mat.zeros D.n 1) ()
let us = Arr.squeeze us and rs = Arr.squeeze rs
let _ = Mat.save_txt (Mat.transpose us) (Dir.in_dir "us_input_steps")
let _ = Mat.save_txt (Mat.transpose rs) (Dir.in_dir "rs_input_steps")



(* --------------------------------------------------------------------------------
   ---    panel G : tuning of variability suppression                           ---
   -------------------------------------------------------------------------------- *)

let contrasts = [ 0.0; 0.5; 1.0 ]

let ff_window = 0.1 (* 100 ms *)
let n_trials = 1000
let burn_in = 1.0
let sample_from = burn_in -. 0.1 (* sample only the last 100 ms *)
let neuron_slice = L [0; 5; 10; 15; 20; 25; 30; 35; 40; 45; 0] 

let prepare v = 
  let v = Arr.transpose ~axis:[| 1; 0; 2 |] v in (* neurons first *)
  let d = Arr.shape v in 
  Arr.reshape v [| d.(0); d.(1) * d.(2) |]

let u0 = Mat.zeros D.n n_trials

(*
let _ = List.iter (fun c ->
    printf "Contrast: %.1f\n%!" c;
    let in_dir s = Dir.in_dir (sprintf "%s_c%.1f" s c) in
    let input = D.(Steady (stim c)) in 
    let us, rs = D.simulate 
        ~verbose:true 
        ~slice:neuron_slice 
        ~duration:burn_in 
        ~input 
        ~u0 () in
    let us' = prepare us and rs' = prepare rs in
    Gc.compact ();
    Mat.(save_txt (mean ~axis:1 us') (in_dir "u_mean"));
    Mat.(save_txt (mean ~axis:1 rs') (in_dir "r_mean"));
    Mat.(save_txt (var ~axis:1 us') (in_dir "u_var"));
    Mat.(save_txt (var ~axis:1 rs') (in_dir "r_var"));
    Mat.(save_txt (D.ff [0;-1] rs) (in_dir "fano_factor"));
  ) contrasts
*)

(* --------------------------------------------------------------------------------
   ---    panel E : variability transient                                       ---
   -------------------------------------------------------------------------------- *)

let n_trials = 2000
let neuron_slice = L [0; 10; 25 ] 

let us, rs = D.simulate 
    ~verbose:true ~duration:1.4
    ~slice:neuron_slice 
    ~sample_from:0.7
    ~input:(transient_input 1.) ~u0:(Mat.zeros D.n n_trials) ()

let u_mean = Arr.(mean ~axis:0 us |> squeeze) |> Mat.transpose
let r_mean = Arr.(mean ~axis:0 rs |> squeeze) |> Mat.transpose
let u_var = Arr.(var ~axis:0 us |> squeeze) |> Mat.transpose
let r_var = Arr.(var ~axis:0 rs |> squeeze) |> Mat.transpose
let ff = D.ff_timecourse ~window:0.1 rs 

let _ =
  Mat.save_txt u_mean (Dir.in_dir "transient_u_mean");
  Mat.save_txt r_mean (Dir.in_dir "transient_r_mean");
  Mat.save_txt u_var (Dir.in_dir "transient_u_var");
  Mat.save_txt r_var (Dir.in_dir "transient_r_var");
  Mat.save_txt ff (Dir.in_dir "transient_fano_factor")





