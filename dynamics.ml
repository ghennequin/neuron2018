open Printf
open Owl
module P = Pervasives

(* ----------------------------------------------------------------------------
   ---   SSN dynamics / sample-based variability estimation                 ---
   ---------------------------------------------------------------------------- *)


module type Prms = sig
  val dt: float
  val sampling_dt: float
  val w: Mat.mat
  val n_e: int
  val n_i: int
  val tau_e: float
  val tau_i: float
  val tau_noise: float
  val sigma_noise: Mat.mat (* noise covariance *)
  val nonlin_k: float (* k of the threshold powerlaw nonlinearity *)
  val nonlin_n: float (* exponent of the threshold powerlaw nonlinearity *)
end


module Make (X: Prms) = struct

  include X

  let n, _ = Mat.shape w
  let _ = assert (n = n_e + n_i)
  let taus = Mat.init n 1 (fun i -> if i<n_e then tau_e else tau_i)

  let prepare_noise n_trials =
    let ell = Linalg.D.chol ~upper:false sigma_noise in
    let eta = ref Mat.(ell *@ gaussian n n_trials) in
    let decay = 1. -. dt /. tau_noise in
    let factor = sqrt (2. *. dt /. tau_noise) in
    fun () ->
      eta := Mat.((decay $* !eta) + (factor $* (ell *@ gaussian n n_trials)));
      !eta

  type input = 
    | Steady of Mat.mat
    | Time_varying of (float -> Mat.mat)

  (* number of columns in initial condition u0 sets the number of trials executed in parallel *)
  let simulate ?(verbose=false) ?slice ?sample_from ~duration ~input ~u0 () =
    let _, n_trials = Mat.shape u0 in
    let n_time_bins = Maths.round (duration /. dt) |> int_of_float in
    let sample_every = Maths.round (sampling_dt /. dt) |> int_of_float in
    let n_samples = match sample_from with
      | None -> n_time_bins / sample_every 
      | Some s -> (Maths.round ((duration -. s) /. dt) |> int_of_float) / sample_every in
    let noise = prepare_noise n_trials in
    let n_slice = match slice with
      | Some s -> fst Mat.(shape (get_fancy [s; R[]] (empty n 1)))
      | None -> n in
    let us = Arr.zeros [| n_samples; n_slice; n_trials |] in
    let rs = Arr.zeros [| n_samples; n_slice; n_trials |] in
    let u = ref u0 in
    let decay = Mat.(1. $- (dt $/ taus)) in
    let dt_over_taus = Mat.(dt $/ taus) in
    for t=0 to n_time_bins-1 do
      let time = dt *. float t in
      if t mod 20 = 0 then Gc.full_major ();
      if t mod 10 = 0 then begin
        if verbose then printf "\r[dynamics] %.3f / %.3f%!" time duration;
      end;
      let r = Mat.(nonlin_k $* (pow_scalar !u nonlin_n)) in
      (* sample *)
      if ((t mod sample_every = 0) 
          && (match sample_from with None -> true | Some s -> (time > s)))
      then begin
        let bin = match sample_from with 
          | None -> t / sample_every 
          | Some s -> (t - int_of_float (Maths.round (s /. dt))) / sample_every in
        Arr.(copy_to 
               (match slice with Some s -> get_fancy [s; R[]] !u | None -> !u) 
               (slice_left us [| bin |]));
        Arr.(copy_to 
               (match slice with Some s -> get_fancy [s; R[]]  r | None ->  r) 
               (slice_left rs [| bin |]));
      end;
      (* update *)
      let h = match input with
        | Steady z -> z
        | Time_varying f -> f time in
      let input = Mat.(h + (w *@ r) + (noise ())) in
      u := Mat.((decay * !u) + (dt_over_taus * input));
    done;
    if verbose then print_newline ();
    let us = Arr.transpose ~axis:[| 2; 1; 0 |] us in
    Gc.full_major ();
    let rs = Arr.transpose ~axis:[| 2; 1; 0 |] rs in
    Gc.full_major ();
    us, rs

  (* spike counts statistics *)

  (* Fano factor; rs: firing rates, trials x neurons x time *)
  let ff time_slice rs =
    let rs = Arr.(rs.${[[]; []; time_slice]}) in
    let kappas = Arr.(sampling_dt $* sum ~axis:2 rs) |> Arr.squeeze (* trials; neurons *) in
    Mat.((1.0 $+ ((var ~axis:0 kappas) / (mean ~axis:0 kappas))) |> transpose)

  (* Fano factor timecourse; rs: firing rates, trials x neurons x time *)
  let ff_timecourse ~window rs =
    let window = Maths.(round (window /. sampling_dt) |> int_of_float) in
    let d = Arr.shape rs in
    let n_time_bins = d.(2) in 
    Array.init (n_time_bins - window) (fun t -> ff [ t; t+window ] rs)
    |> Mat.concatenate ~axis:1 
    |> Mat.transpose

end






