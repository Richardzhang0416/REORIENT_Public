module B_structs
using  O_structs
export D__type    # structure of problem data

mutable struct D__type
   ms::ms_type    # data of master problem sets
   mp::mp_type    # data of master problem parameters
   ps::ps_type    # data of subproblem sets
   pp::pp_type    # data of subproblem parameters
  unc::u_type     # data of uncertain parameters
 algm::Int64      # algorithm used
   γs::Float64    # preset parameter for level method, 0 gives classical method
    J::Int64      # maximum number of iteration
    δ::Float64    # convergence tolerance
end

end
