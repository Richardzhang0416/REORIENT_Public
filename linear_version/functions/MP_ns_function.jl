include("load_stuff.jl")
unc,ms,mp,ps,pp=load_data()
function nxt(t,TB,TE)
    if t< TE
        return t+1
    end
    if t>= TE
        return TB
    end
end
function MP!(ms::ms_type,mp::mp_type,unc::u_type,ps::ps_type,pp::pp_type)
    m=Model(Gurobi.Optimizer)

    # */ ------------------ investment variables --------------------------------/* #
    @variable(m, f>=0) # investment-only cost
    @variable(m, x0ph[ms.PH,ms.ZH,ms.I0]>=0)   # newly installed units capacity on energy hub
    @variable(m, x0pp[ms.PP,ms.ZP,ms.I0]>=0)   # newly installed units capacity on production platform
    @variable(m, x0po[ms.PO,ms.ZO,ms.I0]>=0)   # newly installed units capacity onshore
    @variable(m, x0l[ms.L,ms.I0]>=0)    # newly installed line capacity
    @variable(m, x[unc.tx,ms.I]>=0)    # rhs paramters for operational subproblem (e.g., accumulated capacity)
    @variable(m, β[ms.I]>=0)    # operational cost of node i

    # */ ----------------- operational variables --------------------------- /* #
    @variable(m, ϕ[unc.tc,ms.I]>=0)    # operational cost dependent on uncertain cost cᵢ (now is certain)
    @variable(m, c0[ms.I]>=0)  # operational cost independent of certain costs c

    # */ ---------------- electricity variables ---------------------------- /* #
    @variable(m, pD[ps.ZP,ps.sH,ms.I]>=0)            # total power demand on platform z in period t
    @variable(m, pG[ps.G,ps.ZP,ps.sH,ms.I]>=0)       # power output of gas turbine g in period t
    @variable(m, pResG[ps.G,ps.ZP,ps.sH,ms.I]>=0)    # power reserved for spinning reserve of gas turbine g on platform z in period t
    @variable(m, pSEP[ps.Dsep,ps.ZP,ps.sH,ms.I]>=0)  # power required by separator on platform z in period t
    @variable(m, pRW[ps.RW,ps.ZH,ps.sH,ms.I]>=0)     # wind power output on platform z in period t
    @variable(m, pRS[ps.RS,ps.ZH,ps.sH,ms.I]>=0)     # solar power output on platform z in period t
    @variable(m, pSEI[ps.SE,ps.ZP,ps.sH,ms.I]>=0)    # power injection to electricity store e on platform z in period t
    @variable(m, pSEO[ps.SE,ps.ZP,ps.sH,ms.I]>=0)    # power withdrawn from selectricity store e on platform z in period t
    @variable(m, pResSE[ps.SE,ps.ZP,ps.sH,ms.I]>=0)  # power reserved in electricity store e on platform z in period t
    @variable(m, qSE[ps.SE,ps.ZP,ps.sH,ms.I]>=0)     # energy storage level in electricity store e on platform z in period t
    @variable(m, pGShed[ps.Z,ps.sH,ms.I]>=0)         # generation shed on platform z in period t
    @variable(m, pLShed[ps.ZP,ps.sH,ms.I]>=0)        # load shed on platform z in period t
    @variable(m, pZO[ps.B,ps.ZO,ps.sH,ms.I]>=0)      # power output of onshore bus
    @variable(m, pL[ps.L,ps.sH,ms.I])                # power flow in line l in period t

    # */ ------------------- heat variables -------------------------------- /* #
    @variable(m, pEB[ps.EB,ps.ZP,ps.sH,ms.I]>=0) # heat power output from electric boiler on platform z in period t
    @variable(m, pHGShed[ps.ZP,ps.sH,ms.I]>=0)   # heat generation shed on platform z in period t
    @variable(m, pHLShed[ps.ZP,ps.sH,ms.I]>=0)   # heat load shed on platfom z in period t
    @variable(m, pSEPH[ps.Dsep,ps.ZP,ps.sH,ms.I]>=0) # heat power required by separator

    # */ ------------------ natural gas variables -------------------------- /* #
    @variable(m, pC[ps.C,ps.ZP,ps.sH,ms.I]>=0) # power required by recompressor
    @variable(m, pGI[ps.GI,ps.ZP,ps.sH,ms.I]>=0) # power required by gas injection compressor

    # */ ------------------ water variables -------------------------- /* #
    @variable(m, pO[ps.PO,ps.ZP,ps.sH,ms.I]>=0) # power required by oil pump

    # */ ------------------ oil variables -------------------------- /* #
    @variable(m, pWI[ps.PWI,ps.ZP,ps.sH,ms.I]>=0) # power required by water injection pump
    @variable(m, pWL[ps.PWL,ps.ZP,ps.sH,ms.I]>=0) # power required by water lift pump

    # */ ------------------ hydrogen variables ---------------------------- /* #
    @variable(m, pF[ps.F,ps.ZH,ps.sH,ms.I]>=0)   # power output of fuel cell f on offshore energy hub zh in period t
    @variable(m, pE[ps.E,ps.ZH,ps.sH,ms.I]>=0)   # power input to electrolyser on offshore energy hub zh in period t
    @variable(m, vEF[ps.E,ps.ZH,ps.sH,ms.I]>=0)  # hydrogen produced from electrolyser and supply fuel cell on offshore energy hub zh in period t
    @variable(m, vSHy[ps.SHy,ps.ZH,ps.sH,ms.I]>=0)      # storage level of hydrogen storage on offshore energy hub zh in period t
    @variable(m, vSHyI[ps.SHy,ps.ZH,ps.sH,ms.I]>=0)     # hydrogen injection into storage facility on offshore energy hub zh in period t
    @variable(m, vSHyO[ps.SHy,ps.ZH,ps.sH,ms.I]>=0)     # hydrogen withdrawn from the storage facility on offshore energy hub zh in period t
    @variable(m, x0[unc.tx,ms.I])                # values of rhs parameters in the subproblem (e.g., capacity)

    # */ ------------------ objective function ----------------------------- /* #
    @objective(m, Min, f + sum(mp.π[i]*β[i] for i in ms.I) ) # investment plus operational cost (10⁶£)


    # */ ------------------ master constraints ----------------------------- /* #
    @constraint(m, investment_cost, f>=exp10(-6)*sum(mp.π0[i0]*(sum(mp.ci[mp.map_php[p],i0]*x0ph[p,z,i0] for z in ms.ZH for p in ms.PH) # investment cost
                                                    +sum(mp.ci[mp.map_ppp[p],i0]*x0pp[p,z,i0] for z in ms.ZP for p in ms.PP)
                                                    +sum(mp.ci[mp.map_pop[p],i0]*x0po[p,z,i0] for z in ms.ZO for p in ms.PO)
                                                    +sum(mp.ci[mp.map_ll[l],i0]*x0l[l,i0] for l in ms.L)) for i0 in ms.I0)
                                      +exp10(-6)*sum(mp.π[i]*(sum(mp.cf[mp.map_ppp[p]]*x[(p,z),i] for z in ms.ZP for p in ms.PP)    # fixed O&M cost
                                                    +sum(mp.cf[mp.map_php[p]]*x[(p,z),i] for z in ms.ZH for p in ms.PH)
                                                    +sum(mp.cf[mp.map_pop[p]]*x[(p,z),i] for z in ms.ZO for p in ms.PO)
                                                    +sum(mp.cf[mp.map_ppp[p]]*x[(p,z),i] for z in ms.ZP for p in ms.PP)
                                                    +sum(mp.cf[mp.map_ll[l]]*x[(l,l),i] for l in ms.L)) for i in ms.I))
    @constraint(m, accu_cap_zh[p=ms.PH,z=ms.ZH,i=ms.I], x[(p,z),i]==mp.xhph[mp.map_ph[p],mp.map_zh[z],i]+sum(x0ph[p,z,i0] for i0 in mp.map[i]))
    @constraint(m, accu_cap_zp[p=ms.PP,z=ms.ZP,i=ms.I], x[(p,z),i]==mp.xhpp[mp.map_pp[p],mp.map_zp[z],i]+sum(x0pp[p,z,i0] for i0 in mp.map[i]))
    @constraint(m, accu_cap_zo[p=ms.PO,z=ms.ZO,i=ms.I], x[(p,z),i]==mp.xhpo[mp.map_po[p],mp.map_zo[z],i]+sum(x0po[p,z,i0] for i0 in mp.map[i]))
    @constraint(m, accu_cap_line[l=ms.L,i=ms.I], x[(l,l),i] == mp.xhl[mp.map_l[l],i] + sum(x0l[l,i0] for i0 in mp.map[i]) )
    @constraint(m, max_cap_zh[p=ms.PH,z=ms.ZH,i=ms.I], sum(x0ph[p,z,i0] for i0 in mp.map[i])<=mp.xm0[mp.map_php[p]])
    @constraint(m, max_cap_zp[p=ms.PP,z=ms.ZP,i=ms.I], sum(x0pp[p,z,i0] for i0 in mp.map[i])<=mp.xm0[mp.map_ppp[p]])
    @constraint(m, max_cap_zo[p=ms.PO,z=ms.ZO,i=ms.I], sum(x0po[p,z,i0] for i0 in mp.map[i])<=mp.xm0[mp.map_pop[p]])
    @constraint(m, max_cap_line[l=ms.L,i=ms.I],sum(x0l[l,i0] for i0 in mp.map[i])<=mp.xm0[mp.map_ll[l]])
    for i in ms.I for j in unc.th
        fix(x[(j,"-"),i],unc.h[i,unc.map_h[j]];force=true)   # set rhs paramters to dependent on subproblem i (e.g., CO₂ budget parameter)
    end end

    # */ ---------------- subproblems constraints --------------------------- /* #
    @constraint(m, cost_emission[i=ms.I], β[i] >= c0[i] + unc.c[i,unc.map_c["co2_price"]]*ϕ["co2_price",i]) # compute generation cost (10⁶£) of subproblem i
    @constraint(m, cost_operation[i=ms.I], c0[i] == exp10(-6)*(sum(sum((sum((pp.CG+pp.CFuel)/pp.ηG*pG[g,z,h,i] for g in ps.G)+
                                           pp.CLShed*pLShed[z,h,i]+pp.CHLShed*pHLShed[z,h,i])*pp.Ht for h in ps.sH)*pp.WS for z in ps.ZP)+
                                           sum(sum(sum(pp.CZOE[pp.map_zo[z]]*pZO[b,z,h,i]*pp.EG/pp.ηG for b in ps.B)*pp.Ht for h in ps.sH)*pp.WS for z in ps.ZO)))

    # */ ----------------------- CO_2 emission--------------------------- /* #
    @constraint(m, emisison[i=ms.I], ϕ["co2_price",i] == exp10(-6)*sum(sum(sum(pG[g,z,h,i]*pp.EG/pp.ηG for g in ps.G)*pp.Ht for h in ps.sH)*pp.WS for z in ps.ZP))

    # */ ---------------- electricity constraints ---------------------------------- /* #
    @constraint(m, onshore_bus_cap[b=ps.B, z=ps.ZO, s=ps.S, h=ps.Hs[s],i=ms.I], pZO[b,z,h,i] <= x0[(b,z),i])
    @constraint(m, turbine_cap[g=ps.G, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pG[g,z,h,i]+pResG[g,z,h,i] <= x0[(g,z),i])
    @constraint(m, turbine_ramp_up[g=ps.G, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pG[g,z,h,i]+pResG[g,z,h,i]-pG[g,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]-pResG[g,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]
                <=pp.GGR*x0[(g,z),i])
    @constraint(m, turbine_ramp_down[g=ps.G, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pG[g,z,h,i]+pResG[g,z,h,i]-pG[g,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]-pResG[g,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]
                >=-pp.GGR*x0[(g,z),i])
    @constraint(m, solar[r=ps.RS, z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pRS[r,z,h,i]==pp.RPS[h,pp.map_zh[z]]*x0[(r,z),i])
    @constraint(m, wind[r=ps.RW, z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pRW[r,z,h,i]==pp.RPW[h,pp.map_zh[z]]*x0[(r,z),i])
    @constraint(m, e_store_balance[se=ps.SE, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], qSE[se,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i] == qSE[se,z,h,i]+pp.Ht*(pp.ηSE*pSEI[se,z,h,i]-pSEO[se,z,h,i]))
    @constraint(m, e_store_charge_cap[se=ps.SE, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pSEI[se,z,h,i] <= x0[(se,z),i]*pp.ΗSE)
    @constraint(m, e_store_dicharge_cap[se=ps.SE, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pSEO[se,z,h,i]+pResSE[se,z,h,i] <= x0[(se,z),i]*pp.ΗSE)
    @constraint(m, e_store_energy_cap[se=ps.SE, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], qSE[se,z,h,i] <= x0[(se,z),i])
    @constraint(m, e_store_discharge_energy_relation[se=ps.SE, z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pp.Ht*(pResSE[se,z,h,i]+pSEO[se,z,h,i]) <= qSE[se,z,h,i])
    @constraint(m, KCL_zp[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pG[g,z,h,i] for g in ps.G) +
                    + sum(pp.ηL[pp.map_l[l]]*pp.AE[pp.map_zpz[z],pp.map_l[l]]*pL[l,h,i] for l in ps.L)+ sum(pSEO[se,z,h,i] for se in ps.SE) + pLShed[z,h,i]
                    == pD[z,h,i] + pGShed[z,h,i] + sum(pSEI[se,z,h,i] for se in ps.SE))
    @constraint(m, KCL_zh[z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pRS[r,z,h,i] for r in ps.RS) + sum(pRW[r,z,h,i] for r in ps.RW)
                    + sum(pp.ηL[pp.map_l[l]]*pp.AE[pp.map_zhz[z],pp.map_l[l]]*pL[l,h,i] for l in ps.L) + sum(pF[f,z,h,i] for f in ps.F)
                    == pGShed[z,h,i] + sum(pE[e,z,h,i] for e in ps.E))
    @constraint(m, KCL_zo[z=ps.ZO, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pp.ηL[pp.map_l[l]]*pp.AE[pp.map_zoz[z],pp.map_l[l]]*pL[l,h,i] for l in ps.L) + sum(pZO[b,z,h,i] for b in ps.B)
                    == pGShed[z,h,i])
    @constraint(m, line_cap1[l=ps.L, s=ps.S, h=ps.Hs[s],i=ms.I], pL[l,h,i] <= x0[(l,l),i])
    @constraint(m, line_cap2[l=ps.L, s=ps.S, h=ps.Hs[s],i=ms.I], pL[l,h,i] >= -x0[(l,l),i])
    @constraint(m, reserve[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pp.σRes*pD[z,h,i] <= sum(pResG[g,z,h,i] for g in ps.G) + sum(pResSE[se,z,h,i] for se in ps.SE))
    @constraint(m, electricity_demand[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pD[z,h,i] == sum(pEB[eb,z,h,i] for eb in ps.EB) +
                    sum(pC[c,z,h,i] for c in ps.C) + sum(pO[po,z,h,i] for po in ps.PO) + sum(pGI[gi,z,h,i] for gi in ps.GI) +
                    sum(pWI[pwi,z,h,i] for pwi in ps.PWI) + sum(pWL[pwl,z,h,i] for pwl in ps.PWL))

    # */ ------------------------- heat contraints----------------------------------- /* #
    @constraint(m, E_boiler_cap[eb=ps.EB,z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], pEB[eb,z,h,i] <= x0[(eb,z),i])
    @constraint(m, heat_balance[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pp.ηGh*pG[g,z,h,i] for g in ps.G) + sum(pp.ηEB*pEB[eb,z,h,i] for eb in ps.EB) + pHLShed[z,h,i]
                        == sum(pSEPH[dsep,z,h,i] for dsep in ps.Dsep) + pHGShed[z,h,i])
    @constraint(m, separator[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pSEPH[dsep,z,h,i] for dsep in ps.Dsep) == pp.ρSEP*(pp.VOD[h,pp.map_zp[z]])*x0[("demand_scaling","-"),i])

    # */ ------------------------- natural gas contraints---------------------------- /* #
    @constraint(m, comrepssor[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pC[c,z,h,i] for c in ps.C) == (max((pp.VD[h,pp.map_zp[z]]-pp.VGI[h,pp.map_zp[z]]),0)*0.0051*pp.θNG))
    @constraint(m, inj_compressor[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pGI[gi,z,h,i] for gi in ps.GI) == pp.VGI[h,pp.map_zp[z]]*0.011*pp.θNG)

    # */ ------------------------- oil contraints----------------------------------- /* #
    @constraint(m, power_oil_pump[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pO[po,z,h,i] for po in ps.PO) == pp.κPO*pp.VOD[h,pp.map_zp[z]]*x0[("demand_scaling","-"),i])

    # */ ------------------------- water contraints----------------------------------- /* #
    @constraint(m, power_water_inj[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pWI[pwi,z,h,i] for pwi in ps.PWI) == pp.κPWI*pp.VWI[h,pp.map_zp[z]])
    @constraint(m, power_water_lift[z=ps.ZP, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pWL[pwl,z,h,i] for pwl in ps.PWL) == pp.κPWL*pp.VWL[h,pp.map_zp[z]])

    # */ ------------------------- hydrogen contraints--------------------------------/* #
    @constraint(m, electrolyser_cap[e=ps.E, z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pE[e,z,h,i] <= x0[(e,z),i])
    @constraint(m, fuel_cell_cap[f=ps.F, z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pF[f,z,h,i] <= x0[(f,z),i])
    @constraint(m, hydrogen_strore_cap[shy=ps.SHy,z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], vSHy[shy,z,h,i] <= x0[(shy,z),i])
    @constraint(m, fuel_cell_ramp_up[f=ps.F,z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pF[f,z,h,i]-pF[f,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]<=pp.FR*x0[(f,z),i])
    @constraint(m, fuel_cell_ramp_down[f=ps.F,z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pF[f,z,h,i]-pF[f,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]>=-pp.FR*x0[(f,z),i])
    @constraint(m, hydrogen_store_balance[shy=ps.SHy,z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], vSHy[shy,z,nxt(h,ps.Hs[s][1],ps.Hs[s][end]),i]- vSHy[shy,z,h,i]==
                    vSHyI[shy,z,h,i]-vSHyO[shy,z,h,i])
    @constraint(m, hydrogen_balance1[z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], pp.Ht*sum(pE[e,z,h,i] for e in ps.E) == pp.ηES*sum(vSHyI[shy,z,h,i] for shy in ps.SHy) + sum(vEF[e,z,h,i] for e in ps.E)*pp.ηEF)
    @constraint(m, hydrogen_balance2[z=ps.ZH, s=ps.S, h=ps.Hs[s],i=ms.I], sum(pp.Ht*pF[f,z,h,i]/(pp.ηF*pp.θHY) for f in ps.F) == sum(vEF[e,z,h,i] for e in ps.E) + sum(vSHyO[shy,z,h,i] for shy in ps.SHy))

    # */ ----------------------------- CO2 budget ------------------------------------ /* #
    @constraint(m, co2_budget[i=ms.I], (sum(pG[g,z,h,i]*pp.EG/pp.ηG for g in ps.G for h in ps.sH for z in ps.ZP))*pp.Ht*pp.WS <=
                        pp.ΜE)

    @constraint(m, pass_cap[n=1:length(unc.tx)-2,i=ms.I], x0[unc.tx[n],i]*exp10(-3) == x[unc.tx[n],i])
    @constraint(m, pass_limits[n=unc.th,i=ms.I], x0[(n,"-"),i]   == x[(n,"-"),i])

    return m
end

p=MP!(ms,mp,unc,ps,pp)
optimize!(p)

## output results to directory "plots"
# i=0;investment_cost=[];operational_cost=[];emission=[];co2tax=[];energy_loss=[];location=[];tech=[];line=[];electricity_shed=[];heat_shed=[];co2tax2=[];battery_cap=[];
# wind_cap=[];solar_cap=[];electrolyser_cap=[];fuel_cell_cap=[];hydrogen_store_cap=[];converter_cap=[];co2tax3=[];line_cap=[];line=[];location_EH=[];E_boiler_cap=[];
# co2tax4=[];power_demand=[];heat_demand=[];wind_power=[];gas_power=[];gas_heat=[];onshore_power=[];fuel_cell_power=[];co2tax5=[]; location_ZO=[]; onshore_cap=[];co2budget=[];
# co2budget2=[];co2budget3=[];co2budget4=[];co2budget5=[];el_battery=[];el_boiler=[];el_separator=[];el_inj_comp=[];el_pump_oil=[];el_comp=[];el_wi=[];el_wl=[];el_el=[];el_fc=[];
# el_gs1=[];el_gs2=[];el_gs3=[];el_gs4=[];el_gs5=[];el_gs6=[];el_gs7=[];el_gs8=[];el_gs9=[];el_gs10=[];el_gs11=[];el_gs12=[];el_gs13=[];el_gs14=[];el_gs15=[];
# el_hgs1=[];el_hgs2=[];el_hgs3=[];el_hgs4=[];el_hgs5=[];el_OCGT=[];el_l=[];
# while i<=250
#     i=i+1
#     setfield!(unc,:c,[20.0*(i-1) 20.0*(i-1)])
#     p=MP!(ms,mp,unc,ps,pp)
#     optimize!(p)
#
#         push!(co2tax,20.0*(i-1))
#         push!(investment_cost,value.(p[:f]))
#         push!(operational_cost, value.(p[:β][1]))
#         push!(emission, value.( p[:ϕ]["co2_price",1]))
#         push!(energy_loss, ((sum((1/pp.ηG-1-pp.ηGh)value.(p[:pG]["OCGT",z,h,1])
#         +(1-pp.ηSE)value.(p[:pSEI]["battery",z,h,1])
#         +(1-pp.ηEB)value.(p[:pEB]["E_boiler",z,h,1])
#         +(1-0.9)value.(p[:pSEPH]["separator",z,h,1])
#         +(1-0.8)value.(p[:pGI]["gas_inj_compressor",z,h,1]+p[:pC]["compressor",z,h,1])
#         +(1-0.65)value.(p[:pO]["oil_exp_pump",z,h,1])
#         +(1-0.75)value.(p[:pWI]["water_inj_pump",z,h,1]+p[:pWL]["water_lift_pump",z,h,1])
#         +value.(p[:pHGShed][z,h,1]) for z in ps.ZP for h in ps.sH)
#         +sum(value.(p[:pE]["electrolyser",z,h,1]-pp.θHY*(p[:vEF]["electrolyser",z,h,1]+p[:vSHyI]["hydrogen_store",z,h,1])
#         + (1/pp.ηF-1)*p[:pF]["fuel_cell",z,h,1]) for z in ps.ZH for h in ps.sH))*pp.WS
#         + sum(value.(p[:pGShed][z,h,1] for z in ps.Z for h in ps.sH))*pp.WS
#         + sum(abs((1-pp.ηL[pp.map_l[l]])*value.(p[:pL][l,h,1])) for l in ps.L for h in ps.sH)*pp.WS)/1000) # in GWh
#
#         push!(el_OCGT, (sum((1/pp.ηG-1-pp.ηGh)value.(p[:pG]["OCGT",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_battery, (sum((1-pp.ηSE)value.(p[:pSEI]["battery",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_boiler, (sum((1-pp.ηEB)value.(p[:pEB]["E_boiler",z,h,1]) for z in ps.ZP  for h in ps.sH)*pp.WS)/1000)
#         push!(el_separator, (sum((1-0.9)value.(p[:pSEPH]["separator",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_inj_comp, (sum((1-0.8)value.(p[:pGI]["gas_inj_compressor",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_comp, (sum((1-0.8)value.(p[:pC]["compressor",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_pump_oil, (sum((1-0.65)value.(p[:pO]["oil_exp_pump",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_wi, (sum((1-0.75)value.(p[:pWI]["water_inj_pump",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_wl, (sum((1-0.75)value.(p[:pWL]["water_lift_pump",z,h,1]) for z in ps.ZP for h in ps.sH)*pp.WS)/1000)
#         push!(el_hgs1, sum(value.(p[:pHGShed]["cluster1",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_hgs2, sum(value.(p[:pHGShed]["cluster2",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_hgs3, sum(value.(p[:pHGShed]["cluster3",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_hgs4, sum(value.(p[:pHGShed]["cluster4",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_hgs5, sum(value.(p[:pHGShed]["cluster5",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs1, sum(value.(p[:pGShed]["cluster1",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs2, sum(value.(p[:pGShed]["cluster2",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs3, sum(value.(p[:pGShed]["cluster3",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs4, sum(value.(p[:pGShed]["cluster4",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs5, sum(value.(p[:pGShed]["cluster5",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs6, sum(value.(p[:pGShed]["EH1",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs7, sum(value.(p[:pGShed]["EH2",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs8, sum(value.(p[:pGShed]["EH3",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs9, sum(value.(p[:pGShed]["EH4",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs10, sum(value.(p[:pGShed]["EH5",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs11, sum(value.(p[:pGShed]["onshore1",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs12, sum(value.(p[:pGShed]["onshore2",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs13, sum(value.(p[:pGShed]["onshore3",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs14, sum(value.(p[:pGShed]["onshore4",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_gs15, sum(value.(p[:pGShed]["onshore5",h,1]) for h in ps.sH)*pp.WS/1000)
#         push!(el_el,sum(value.(p[:pE]["electrolyser",z,h,1]-pp.θHY*(p[:vEF]["electrolyser",z,h,1]+p[:vSHyI]["hydrogen_store",z,h,1])) for h in ps.sH for z in ps.ZH)*pp.WS/1000)
#         push!(el_fc, sum((1/pp.ηF-1)*value.(p[:pF]["fuel_cell",z,h,1]) for z in ps.ZH for h in ps.sH)*pp.WS/1000)
#         push!(el_l, sum(abs((1-pp.ηL[pp.map_l[l]])*value.(p[:pL][l,h,1])) for l in ps.L for h in ps.sH)*pp.WS/1000)
#
#
#         for z in ms.ZP
#             push!(co2tax2,20.0*(i-1))
#             push!(location,z)
#             push!(electricity_shed,sum(value.(p[:pLShed][z,h,1] for h in ps.sH)))
#             push!(heat_shed,sum(value.(p[:pHLShed][z,h,1] for h in ps.sH))/2880)
#             push!(power_demand, sum(value.(p[:pD][z,h,1] for h in ps.sH))/2880)
#             push!(heat_demand, sum(value.(p[:pSEPH][dsep,z,h,1] for dsep in ps.Dsep for h in ps.sH))/2880)
#             push!(gas_power, sum(value.(p[:pG][g,z,h,1] for g in ps.G for h in ps.sH))/2880)
#             push!(gas_heat, sum(value.(pp.ηGh*p[:pG][g,z,h,1] for g in ps.G for h in ps.sH))/2880)
#             push!(battery_cap,value.(p[:x0pp]["battery",z,1])*exp10(3))
#             push!(E_boiler_cap,value.(p[:x0pp]["E_boiler",z,1])*exp10(3))
#
#         end
#
#         for z in ms.ZH
#             push!(co2tax4,20.0*(i-1))
#             push!(location_EH,z)
#             push!(wind_cap,value.(p[:x0ph]["wind",z,1])*exp10(3))
#             push!(solar_cap,value.(p[:x0ph]["solar",z,1])*exp10(3))
#             push!(electrolyser_cap,value.(p[:x0ph]["electrolyser",z,1])*exp10(3))
#             push!(fuel_cell_cap,value.(p[:x0ph]["fuel_cell",z,1])*exp10(3))
#             push!(hydrogen_store_cap,value.(p[:x0ph]["hydrogen_store",z,1])) # in tonne
#             push!(converter_cap,value.(p[:x0ph]["converter",z,1])*exp10(3))
#             push!(wind_power, sum(value.(p[:pRW][w,z,h,1] for w in ps.RW for h in ps.sH))/2880)
#             push!(fuel_cell_power, sum(value.(p[:pF][f,z,h,1] for f in ps.F for h in ps.sH))/2880)
#
#         end
#
#         for z in ms.ZO
#             push!(co2tax5,20.0*(i-1))
#             push!(location_ZO,z)
#             push!(onshore_power,sum(value.(p[:pZO][b,z,h,1] for b in ps.B for h in ps.sH))/2880)
#             push!(onshore_cap, maximum(value.(p[:pZO][ps.B,z,ps.sH,1])))
#         end
#
#         for l in ps.L
#             push!(co2tax3,20.0*(i-1))
#             push!(line_cap,value.(p[:x0l][l,1])*exp10(3))
#             push!(line,l)
#         end
# end
#
# cost_D=DataFrame(co2tax=co2tax,investment_cost=investment_cost, operational_cost=operational_cost)
# emission_D=DataFrame(co2tax=co2tax, emission=emission, energy_loss=energy_loss, el_battery=el_battery, el_boiler=el_boiler, el_separator=el_separator, el_inj_comp=el_inj_comp,
# el_pump_oil=el_pump_oil,el_comp=el_comp, el_wi=el_wi,el_wl=el_wl,el_el=el_el,el_fc=el_fc,
# el_gs1=el_gs1,el_gs2=el_gs2,el_gs3=el_gs3,el_gs4=el_gs4,el_gs5=el_gs5,el_gs6=el_gs6,el_gs7=el_gs7,el_gs8=el_gs8,el_gs9=el_gs9,el_gs10=el_gs10,el_gs11=el_gs11,el_gs12=el_gs12,
# el_gs13=el_gs13,el_gs14=el_gs14,el_gs15=el_gs15,
# el_hgs1=el_hgs1,el_hgs2=el_hgs2,el_hgs3=el_hgs3,el_hgs4=el_hgs4,el_hgs5=el_hgs5,el_OCGT=el_OCGT,el_l=el_l)
# shed_D=DataFrame(co2tax=co2tax2,location=location,electricity_shed=electricity_shed,heat_shed=heat_shed)
# tech_D=DataFrame(co2tax=co2tax2,location=location_EH,wind=wind_cap,solar=solar_cap,electrolyser=electrolyser_cap,fuel_cell=fuel_cell_cap,hydrogen=hydrogen_store_cap,converter=converter_cap,battery=battery_cap,E_boiler=E_boiler_cap)
# line_D=hcat(DataFrame(line=line),DataFrame(co2tax=co2tax3),DataFrame(line_cap=line_cap))
# power_demand_D=DataFrame(co2tax=co2tax2,power_demand=power_demand,location=location)
# heat_demand_D=DataFrame(co2tax=co2tax2,heat_demand=heat_demand,location=location)
# gas_power_D=DataFrame(co2tax=co2tax2,OCGT=gas_power,location=location)
# power_supply_D=DataFrame(co2tax=co2tax4,wind_power=wind_power,fuel_cell_power=fuel_cell_power,location=location_EH)
# power_onshore_D=DataFrame(co2tax=co2tax5,onshore_power=onshore_power,location=location_ZO,onshore_cap=onshore_cap)
#
# XLSX.writetable(pwd()*"/plots/co2tax_sensitivity_5000.xlsx",overwrite=true,
# cost=(collect(DataFrames.eachcol(cost_D)), DataFrames.names(cost_D)),
# emission=(collect(DataFrames.eachcol(emission_D)), DataFrames.names(emission_D)),
# shed=(collect(DataFrames.eachcol(shed_D)), DataFrames.names(shed_D)),
# tech=(collect(DataFrames.eachcol(tech_D)), DataFrames.names(tech_D)),
# line=(collect(DataFrames.eachcol(line_D)), DataFrames.names(line_D)),
# power_demand=(collect(DataFrames.eachcol(power_demand_D)), DataFrames.names(power_demand_D)),
# heat_demand=(collect(DataFrames.eachcol(heat_demand_D)), DataFrames.names(heat_demand_D)),
# gas_power=(collect(DataFrames.eachcol(gas_power_D)), DataFrames.names(gas_power_D)),
# power_supply=(collect(DataFrames.eachcol(power_supply_D)), DataFrames.names(power_supply_D)),
# power_onshore=(collect(DataFrames.eachcol(power_onshore_D)), DataFrames.names(power_onshore_D)));
