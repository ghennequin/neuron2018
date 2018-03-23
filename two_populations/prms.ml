open Owl

module Default = struct
  include Generic_prms.Default
  let n_e = 1
  let n_i = 1
  let w = Generic_prms.(weight_matrix default_weights)
  let sigma_noise = Generic_prms.sigma_noise 
      ~tau_e ~tau_i ~tau_noise
      ~sigma_e:0.2 ~sigma_i:0.1 ~rho:0.
end

