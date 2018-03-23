open Owl

let deg x = x *. Const.pi /. 180.

type stim_prms = { width: float; 
                   position: float; 
                   contrast: float;
                   baseline: float }

(* direction preference of cell i *)
let pref_dir ~m i =
  assert (i<m); 
  2. *. Const.pi *. float i /. float m

(* Von Mises fall off with angular distance *)
let von_mises ~ell =
  let kappa = 1. /. Maths.sqr ell in
  fun theta1 theta2 -> Maths.(exp (kappa *. (cos (theta1 -. theta2) -. 1.)))

(* builds the ring weight matrix *)
let weight_matrix ~m ~ell prms =
  let theta = pref_dir ~m in
  let falloff = von_mises ~ell in
  let base = 
    let b = Mat.init_2d m m (fun i j -> falloff (theta i) (theta j)) in
    let z = 1. /. Mat.(mean' (sum ~axis:1 b)) in
    Mat.(z $* b) in
  Mat.kron Generic_prms.(weight_matrix prms) base

let stimulus ~m ~ell ~location ~ie_ratio =
  let theta = pref_dir ~m in
  let falloff = von_mises ~ell in
  let base = Mat.init m 1 (fun k -> falloff (theta k) location) in
  Mat.concat_vertical base Mat.(ie_ratio $* base)

let sigma_noise ~tau_e ~tau_i ~tau_noise ~sigma_e ~sigma_i ~rho ~m ~ell =
  let theta = pref_dir ~m in
  let falloff = von_mises ~ell in
  let base = Mat.init_2d m m (fun i j -> falloff (theta i) (theta j)) in
  let twod = Generic_prms.sigma_noise ~tau_e ~tau_i ~tau_noise ~sigma_e ~sigma_i ~rho in
  Mat.add_diag (Mat.kron twod base) 1E-5

module Default = struct
  include Generic_prms.Default
  let m = 50
  let n_e = m
  let n_i = m
  let w = weight_matrix ~m ~ell:(deg 45.) Generic_prms.default_weights
  let sigma_noise = sigma_noise ~tau_e ~tau_i ~tau_noise 
      ~sigma_e:1. ~sigma_i:0.5 ~rho:1.
      ~m ~ell:(deg 60.)
  let stimulus = stimulus ~m ~ell:(deg 60.) ~location:(deg 180.) ~ie_ratio:1.
end

