module O_structs

export u_type   # structure of uncertainty parameters
export ms_type  # structure of master problem sets
export mp_type  # structure of master problem parameters
export ps_type  # structure of subproblem sets
export pp_type  # structure of subproblem parameters

################## structure of uncertainty parameters ############################
mutable struct u_type
    ni::Int64                           # number of subproblems
    tx::Vector{Tuple{String,String}}    # tuple indeices of vector x (technologies, demand and CO2 scaling)
    th::Array{String,1}                 # sets of right hand side vector h (demand and CO2 scaling)
    tc::Array{String,1}                 # sets of cost vector c (CO2 tax)
     h::Array{Float64,2}                # array of right hand side vector h
     c::Array{Float64,2}                # array of cost vector c
 map_h::Dict{String,Int64}              # map th to h
 map_c::Dict{String,Int64}              # map tc to c
end

################## structure of master problem sets ############################
mutable struct ms_type
     P::Array{String,1}       # set of all devices
    PP::Array{String,1}       # set of devices on platform
    PH::Array{String,1}       # set of devices on energy hub
    PO::Array{String,1}       # set of devices onshore
     Z::Array{String,1}       # set of all platforms
    ZP::Array{String,1}       # set of platforms
    ZH::Array{String,1}       # set of energy hub
    ZO::Array{String,1}       # set of onshore nodes
     L::Array{String,1}       # set of lines
    I0::UnitRange{Int64}      # set of "investment" nodes
     I::UnitRange{Int64}      # set of "operational" nodes
end

################## structure of master problem parameters ############################
mutable struct mp_type
    κ::Float64                # years of operational problem
   ##### historical capacity #####
   xhl::Array{Float64,2}        # historical line capacity (GW)
   xhpp::Array{Float64,3}       # historical installed capacity of devices on platform(GW,tonne)
   xhph::Array{Float64,3}       # historical installed capacity of devices on energy hub (GW,tonne)
   xhpo::Array{Float64,3}       # historical installed capacity of devices onshore(GW,tonne)
   ##### unit capacity
   xucap::Array{Float64,1}      # unit capacity of each technology (GW)
   ####### maximum newly inveseted capacity
   xm0::Array{Float64,1}        # maximum newly invested capacity
   ##### cost of devices#####
   ci::Array{Float64,2}        # investment cost
   cf::Array{Float64,1}        # fixed OM cost
  cfi::Array{Float64,1}        # fixed capex
   ##### probabilities #####
   π0::Array{Float64,1}       # probability associated to "investment" node i0 (-)
    π::Array{Float64,1}       # probability associated to "operational" node i (-)
    ##### linking MP and SP #####
    map::Vector{Array{Int64,1}}       # map "operational" node i → "investment" node i0"
map_ppp::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPP... )"
map_php::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPH... )"
map_pop::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
 map_zp::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
 map_zh::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
 map_zo::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
 map_ll::Dict{String,Int64}    # map "(line 𝑙)" -> "(position in cost arrays xl... )"
 map_pp::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPP... )"
 map_ph::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPH... )"
 map_po::Dict{String,Int64}    # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
  map_l::Dict{String,Int64}    # map "(line 𝑙)" -> "(position in cost arrays xl... )"
end


###################### structure of subproblem sets #######################################
mutable struct ps_type

     ######### general ###########
      S::UnitRange{Int64}       # set of time slices
      H::UnitRange{Int64}       # set of all periods
     sH::Array{Int64,1}        # set of all periods in all slices
     Hs::Array{UnitRange{Int64},1} # set of periods per slices
      P::Array{String,1}       # set of all technologies
      Z::Array{String,1}       # set of all technologies
     ZP::Array{String,1}       # set of production platforms
     ZH::Array{String,1}       # set of energy hub platforms
     ZO::Array{String,1}       # set of onshore buses


     ######### electricity ###########
      L::Array{String,1}       # set of transmission lines
      G::Array{String,1}       # set of gas turbines
     RW::Array{String,1}       # set of wind
     RS::Array{String,1}       # set of solar
     SE::Array{String,1}       # set of electricity stores
     EB::Array{String,1}       # set of electricity boiler
   Dcov::Array{String,1}       # set of AC/DC converter
      B::Array{String,1}       # set of onshore buses

      ######### natural gas ###########
      W::Array{String,1}       # set of gas wells
      C::Array{String,1}       # set of gas compressor
     GB::Array{String,1}       # set of gas boilers
     GI::Array{String,1}       # set of gas injection pump
   Dsep::Array{String,1}       # set of separator devices

     ######### oil ###########
     PO::Array{String,1}       # set of oil pump

     ######### water ###########
     PWI::Array{String,1}       # set of water injection pump
     PWL::Array{String,1}       # set of wter lift pump

      ######### hydrogen ###########
      E::Array{String,1}      # set of electrolyser
      F::Array{String,1}      # set of fuel cell
    SHy::Array{String,1}      # set of hydrogen storage
end

################## structure of subproblem parameters ############################
mutable struct pp_type

      ######### general ###########
      Ht::Float64                # number of hours in one period (h)
      WS::Float64                # slices weight (-)

      ######### electricity ###########
     GGR::Float64                # ramping limit of gas turbine (MW/MW)
     RPW::Array{Float64,2}       # wind generation profile (MW/MW)
     RPS::Array{Float64,2}       # solar generation profile (MW/MW)
      ηG::Float64                # efficiency of gas turbine
      ηL::Array{Float64,1}       # efficiency of line
      EG::Float64                # CO2 emissions per MWh power generated (tonne/MWh)
     ηSE::Float64                # efficiency of electricity store
     ΗSE::Float64                # ratio energy to power of electricity store
      AE::Array{Float64,2}       # bus-line incidence matrix
    κSEP::Float64                # electric demand as fraction of amount of gas passed throguh separator
   CFuel::Float64                # fuel cost of gas trubine (kr/MWh)
      CG::Float64                # operational cost of generating 1MW from gas turbine (kr/MW)
  CLShed::Float64                # electricity load shed cost (kr/MW)
  CGShed::Float64                # generation shed cost (kr/MW)
     CSE::Float64                # cost of power into electricity store (kr/MW)
    σRes::Float64                # spinning reserce factor
     CZO::Array{Float64,2}       # electricity price onshore
    CZOP::Array{Float64,1}       # abatement cost using electricity matrix
    CZOE::Array{Float64,1}       # abatement cost using emission matrix

      ######### heat ###########
     ηGh::Float64                # heat recovery efficiency of gas turbine
     ηFh::Float64                # heat recovery efficiency of fuel cell
    ρSEP::Float64                # heat demand as fraction of amount of gas separated in separator
     ηEB::Float64                # efficiency of gas boiler
 CHLShed::Float64                # heat load shedding

      ######### natural gas ###########
      VD::Array{Float64,2}        # natural gas production level (kg)
     VGI::Array{Float64,2}        # natural gas production level (kg)
 CNGShed::Float64                 # natural gas shed cost (kr/kg)
      CC::Float64                 # operational cost of compressor (kr/kg)
    CSEP::Float64                 # operational cost of separator (kr/kg)
      CW::Float64                 # operational cost of well production (kr/kg)
     CGB::Float64                 # operational cost of natural gas boiler (kr/MW)
     ηGB::Float64                 # efficiency of natural gas boiler
      γC::Float64                 # compression ratio of compressor
       α::Float64                 # polytropic exponent of gas compressor
      ηC::Float64                 # efficiency of compressor
     θNG::Float64                 # natural gas energy content (MWh/kg)

     ######### oil ###########
     VOD::Array{Float64,2}        # oil exporting level (kg)
     κPO::Float64                 # electric demand as fraction of amount of oil passed throguh oil exporting pump

     ######### water ###########
     VWB::Array{Float64,2}        # water from bore (kg)
     VWI::Array{Float64,2}        # water injection level (kg)
     VWL::Array{Float64,2}        # water lift level (kg)
     κPWI::Float64                # electric demand as fraction of amount of water injected by water injection pump
     κPWL::Float64                # electric demand as fraction of amount of water lifted by water pump

      ######### hydrogen ###########
      FR::Float64                 # ramping limit of fuel cell (MW/MW)
    CSHy::Float64                 # operational cost of storing hydrogen in storage facility (kr/kg)
      CE::Float64                 # operational cost of turning electricity into hydrogen in electrolyser (kr/MWh)
      CF::Float64                 # operational cost of turning hydrogen to electricity in fuel cell (kr/MWh)
     ηES::Float64                 # conversion factor of electrolyser to inject gas to the storage facility
     ηEF::Float64                 # conversion factor of electrolyser to inject gas directly to fuel cell
      ηF::Float64                 # efficiency of fuel cell
      σE::Float64                 # minimum allowed working power of electrolyser, fraction of rated capacity
    σSHy::Float64                 # cushion gas, fraction of rated storage capacity
     θHY::Float64                 # hydrogen energy content (MWh/kg)

      ######### other ###########
       μD::Float64                # scaling demand (-)
       μE::Float64                # scaling CO2 yearly limit (tCO₂)
       ΜE::Float64                # initial CO2 limit yearly
     CCO2::Float64                # CO2 tax (kr/tonne)

     ####### map
   map_zp::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
   map_zh::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
   map_zo::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
  map_zpz::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
  map_zhz::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
  map_zoz::Dict{String,Int64}     # map "(tech 𝑝)" -> "(position in cost arrays xhPO... )"
    map_l::Dict{String,Int64}     # map "(line 𝑙)" -> "(position in cost arrays xl... )"

end

end
