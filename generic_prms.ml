open Owl

module Default = struct
  let dt = 1E-4
  let sampling_dt = 1E-3
  let tau_e = 20E-3
  let tau_i = 10E-3
  let tau_noise = 50E-3
  let nonlin_k = 0.3
  let nonlin_n = 2.
end

type t_weights = { w_ee: float; w_ie: float; w_ei: float; w_ii: float }
let default_weights = { w_ee = 1.25; w_ie = 1.2; w_ei = -0.65; w_ii = -0.5 }

(* weight matrix *)
let weight_matrix prms = Mat.of_arrays
    [| [| prms.w_ee; prms.w_ei |]; [| prms.w_ie; prms.w_ii |] |]

(* sigma_{e/i}: desired u-std in the unconnected net;
   rho: desired input correlation *)
let sigma_noise ~tau_e ~tau_i ~tau_noise ~sigma_e ~sigma_i ~rho =
  let sigma_e = sigma_e *. sqrt (1. +. tau_e /. tau_noise) in
  let sigma_i = sigma_i*. sqrt (1. +. tau_i /. tau_noise) in
  let s2e = Maths.sqr sigma_e in
  let s2i = Maths.sqr sigma_i in
  let c = rho *. sigma_e *. sigma_i in 
  Mat.of_arrays [| [| s2e; c |]; [| c; s2i |] |]

