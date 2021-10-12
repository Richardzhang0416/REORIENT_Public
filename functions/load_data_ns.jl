function load_data()
    path = pwd()
    XFA = XLSX.readxlsx(path*"/data/algorithm_data_ns.xlsx")
    XFI = XLSX.readxlsx(path*"/data/investment_data_ns.xlsx")
    XFO_Electricity = XLSX.readxlsx(path*"/data/operation_electricity_data_ns.xlsx")
    XFO_heat = XLSX.readxlsx(path*"/data/operation_heat_data_ns.xlsx")
    XFO_NG = XLSX.readxlsx(path*"/data/operation_NG_data_ns.xlsx")
    XFO_W = XLSX.readxlsx(path*"/data/operation_water_data_ns.xlsx")
    XFO_O = XLSX.readxlsx(path*"/data/operation_oil_data_ns.xlsx")
    XFO_hydrogen = XLSX.readxlsx(path*"/data/operation_hydrogen_data_ns.xlsx")
    XFU = XLSX.readxlsx(path*"/data/uncertainty_data_ns.xlsx")

    unc= u_type(0,Vector{Tuple{String,String}}(),String[],String[],zeros(0,0),zeros(0,0),Dict{String,Int64}(),Dict{String,Int64}())
    ms = ms_type(String[],String[],String[],String[],String[],String[],String[],String[],String[],0:0,0:0)
    mp = mp_type(0,
                zeros(0,0),zeros(0,0,0),zeros(0,0,0),zeros(0,0,0),
                zeros(0),
                zeros(0),
                zeros(0,0),zeros(0),zeros(0),
                zeros(0),zeros(0),
                Vector{Array{Int64,1}}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),
                Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),
                Dict{String,Int64}(),Dict{String,Int64}())
    ps = ps_type(0:0,0:0,Array{Int64,1}(),Array{UnitRange{Int64},1}(),String[],String[],String[],String[],String[],
                String[],String[],String[],String[],String[],String[],String[],String[],
                String[],String[],String[],String[],String[],
                String[],
                String[],String[],
                String[],String[],String[])
    pp = pp_type(0.,0.,
                0.,zeros(0,0),zeros(0,0),0.,zeros(0),0.,0.,0.,zeros(0,0),0.,0.,0.,0.,0.,0.,0.,zeros(0,0),zeros(0),zeros(0),
                0.,0.,0.,0.,0.,
                zeros(0,0),zeros(0,0),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                zeros(0,0), 0.,
                zeros(0,0),zeros(0,0),zeros(0,0),0.,0.,
                0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                0.,0.,0.,0.,
                Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}(),Dict{String,Int64}())


    # */ ------------------------- structure ----------------------------------- /* #
    sheet = XFU["structure"]

    ndict = Dict(:t0 =>"present",
                 :t5 =>"5 years")

    tarray = [:t0,:t5]

    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    t0__col = findfirst(col_names .== ndict[:t0 ])[2]
    t5__col = findfirst(col_names .== ndict[:t5 ])[2]
    t0__value = convert(Array{Int64,1},sheet[XLSX.encode_column_number(t0__col)*"2"*":"*XLSX.encode_column_number(t0__col)*"$stpr"][:,1])
    t5__value = convert(Array{Int64,1},sheet[XLSX.encode_column_number(t5__col)*"2"*":"*XLSX.encode_column_number(t5__col)*"$stpr"][:,1])
    nω1 = length(unique(t0__value))
    nω2 = length(unique(t5__value))

    d_dict = Dict{Symbol,Array{Tuple{Int64,UnitRange{Int64},Int64},1}}()
    for sym in keys(ndict)
        tc = findfirst(col_names .== ndict[sym])[2]
        tv = convert(Array{Int64,1},sheet[XLSX.encode_column_number(tc)*"2"*":"*XLSX.encode_column_number(tc)*"$stpr"][:,1])
        ay = [(tc,1:length(tv),tc)]
        d_dict[sym] = ay
    end
    d_array = [(1,2)]
    s_dict = Dict{Symbol,Array{Tuple{Int64,Int64,Int64},1}}()
    for sym in keys(ndict)
        tc = findfirst(col_names .== ndict[sym])[2]
        tv = convert(Array{Int64,1},sheet[XLSX.encode_column_number(tc)*"2"*":"*XLSX.encode_column_number(tc)*"$stpr"][:,1])
        tn = unique(tv)
        ay = Array{Tuple{Int64,Int64,Int64},1}(undef,length(tn))
        for i in 1:length(tn)
            tr = findfirst(tv .== tn[i])
            ay[i] = (tn[i],tr,tc)
        end
        s_dict[sym] = ay
    end
    s_array = Array{Tuple{Int64,Int64},1}(undef,length(t0__value))
    for i in 1:length(t0__value)
        s_array[i] = (t0__value[i],t5__value[i])
    end


    ################ parameter ##################
    sheet = XFU["parameters"]
    ndict = Dict(:shn =>"sheet_name",
                :rc  =>"location",
                :type=>"type")
    udict = Dict(:μD  =>"demand_scaling",
                :μE =>"co2_budget_scaling",
                :CCO2 =>"co2_price")
    uarray = [:μD,:μE,:CCO2]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    row_names = sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"]
    l_col = findfirst(col_names .== ndict[:rc  ])[2]
    t_col = findfirst(col_names .== ndict[:type])[2]
    l_ay = convert(Array{Int64,1},sheet[XLSX.encode_column_number(l_col)*"3"*":"*XLSX.encode_column_number(l_col)*"$stpr"][:,1])
    t_ay = convert(Array{Int64,1},sheet[XLSX.encode_column_number(t_col)*"3"*":"*XLSX.encode_column_number(t_col)*"$stpr"][:,1])
    NI1 = nω1^sum(t_ay)
    NI2 = nω2^sum(t_ay)
    NI0 = NI1
    NI = NI2
    th = String[];tc =String[];
    for i in 1:length(l_ay)
        if l_ay[i]==0
            push!(th,String(row_names[i]))
        elseif l_ay[i]==1
            push!(tc,String(row_names[i]))
    end
    end
    nh = length(l_ay)-sum(l_ay)
    nc = sum(l_ay)


    scol = findfirst(col_names .== ndict[:shn])[2]
    svalue = convert(Array{String,1},sheet[XLSX.encode_column_number(scol)*"3"*":"*XLSX.encode_column_number(scol)*"$stpr"][:,1])
    ddict = Dict{Symbol,Dict{Int64,Float64}}()
    adict = Dict{Symbol,Array{Tuple{Int64,Int64},1}}()

    for sym in uarray
        row = findfirst(row_names .== udict[sym])[1]
        sname = svalue[row]
        sheet = XFU[sname]
        stt,stp = sheet.dimension.start,sheet.dimension.stop
        sttr,sttc = stt.row_number,stt.column_number
        stpr,stpc = stp.row_number,stp.column_number
        values = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(stpc)*"$stpr"])
        vdict = Dict{Int64,Float64}()
        if t_ay[row] == 0
        for tm in tarray
            dic_tm = d_dict[tm]
            for j in 1:length(dic_tm)
                v = values[dic_tm[j][2]]
                vdict[dic_tm[j][1]] = sum(v)/length(v)
            end
        end
        varray = d_array
        else
        for tm in tarray
            dic_tm = s_dict[tm]
            for j in 1:length(dic_tm)
                vdict[dic_tm[j][1]] = values[dic_tm[j][2]]
            end
        end
        varray = s_array
        end
        ddict[sym] = vdict
        adict[sym] = varray
    end

    glb_array = Array{typeof(ntuple(x->0,length(uarray))),2}(undef,0,2)
    for i1 in adict[uarray[1]] for i2 in adict[uarray[2]] for i3 in adict[uarray[3]]
        ay = Array{typeof(ntuple(x->0,length(uarray))),2}(undef,1,2)
        for j in 1:2 ay[j] = (i1[j],i2[j],i3[j]) end
        glb_array = vcat(glb_array,ay)
    end end end
    num_dict = Dict{Int64,typeof(glb_array[1,1])}()
    dict_num = Dict{typeof(glb_array[1,1]),Int64}()
    i_n = 1
    for j in 1:2
        for el in unique(glb_array[:,j])
            num_dict[i_n] = el
            dict_num[el] = i_n
            i_n += 1
        end
    end
    nr = i_n-2
    H,C = zeros(nr,nh),zeros(nr,nc)
    for i in 1:nr
        ih,ic=1,1
        dc = num_dict[i+1]
        for j in 1:length(dc)
            if l_ay[j] == 0
                H[i,ih] = ddict[uarray[j]][dc[j]]
                ih += 1
            end
            if l_ay[j] == 1
                C[i,ic] = ddict[uarray[j]][dc[j]]
                ic += 1
            end
        end
    end

    itoi0 = Vector{Array{Int64,1}}()
    for i in 1:NI1
        rw = findfirst(x->x==num_dict[i+NI1],glb_array[:,2])
        push!(itoi0,[dict_num[glb_array[rw,1]]])
    end
    # for i in 1:NI2
    #     rw = findfirst(x->x==num_dict[i+NI1+NI2],glb_array[:,3])
    #     push!(itoi0,[dict_num[glb_array[rw,1]],dict_num[glb_array[rw,2]]])
    # end

    setfield!(unc,:ni,NI)
    setfield!(unc,:th,th)
    setfield!(unc,:tc,tc)
    setfield!(unc,:h,H)
    setfield!(unc,:c,C)

    map_h=Dict{String,Int64}();map_c=Dict{String,Int64}();
    for i in 1:length(th)
        map_h[th[i]]=i
    end
    for i in 1:length(tc)
        map_c[tc[i]]=i
    end
    setfield!(unc,:map_h,map_h)
    setfield!(unc,:map_c,map_c)
    setfield!(ms,:I0,1:NI0)
    setfield!(ms,:I ,1:NI)
    setfield!(mp,:π0,vcat(ones(NI1)./NI1))
    setfield!(mp,:π ,vcat(ones(NI2)./NI2))
    setfield!(mp,:map,itoi0)

    ########## read nodes sets ##########
    sheet = XFI["nodes"]
    ndict = Dict(:node =>"node",
                :type   =>"type")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    fcol = findfirst(col_names .== ndict[:node])[2]
    tcol = findfirst(col_names .== ndict[:type])[2]
    fvalue = convert(Array{String,1},sheet[XLSX.encode_column_number(fcol)*"3"*":"*XLSX.encode_column_number(fcol)*"$stpr"][:,1])
    tvalue = sheet[XLSX.encode_column_number(tcol)*"3"*":"*XLSX.encode_column_number(tcol)*"$stpr"][:,1]
    ZP=String[];ZH=String[];ZO=String[];Z=String[];
    for i in 1:length(tvalue)
        if tvalue[i]==0
            push!(ZP,fvalue[i])
            push!(Z,fvalue[i])
        elseif tvalue[i]==1
            push!(ZH,fvalue[i])
            push!(Z,fvalue[i])
        elseif tvalue[i]==2
            push!(ZO,fvalue[i])
            push!(Z,fvalue[i])
        end
    end
    setfield!(ms,:ZP,ZP)
    setfield!(ms,:ZH,ZH)
    setfield!(ms,:ZO,ZO)
    setfield!(ms,:Z,Z)
    setfield!(ps,:ZP,ZP)
    setfield!(ps,:ZH,ZH)
    setfield!(ps,:ZO,ZO)
    setfield!(ps,:Z,Z)
    map_zp=Dict{String,Int64}();map_zh=Dict{String,Int64}();map_zo=Dict{String,Int64}();
    for i in 1:length(ms.ZP)
        map_zp[ms.ZP[i]]=i
    end
    for i in 1:length(ms.ZH)
        map_zh[ms.ZH[i]]=i
    end
    for i in 1:length(ms.ZO)
        map_zo[ms.ZO[i]]=i
    end
    map_zpz=Dict{String,Int64}();map_zhz=Dict{String,Int64}();map_zoz=Dict{String,Int64}();
    for i in 1:length(ps.Z)
        if tvalue[i]==0
            map_zpz[ps.Z[i]]=i
        elseif tvalue[i]==1
            map_zhz[ps.Z[i]]=i
        elseif tvalue[i]==2
            map_zoz[ps.Z[i]]=i
        end
    end
    setfield!(mp,:map_zp,map_zp)
    setfield!(mp,:map_zh,map_zh)
    setfield!(mp,:map_zo,map_zo)
    setfield!(pp,:map_zp,map_zp)
    setfield!(pp,:map_zh,map_zh)
    setfield!(pp,:map_zo,map_zo)
    setfield!(pp,:map_zpz,map_zpz)
    setfield!(pp,:map_zhz,map_zhz)
    setfield!(pp,:map_zoz,map_zoz)

    ########### read investment ################
    sheet = XFI["investment"]
    ndict = Dict(:device=>"device",
                :ci=>"CaPeX",
                :cf=>"FixOM",
                :lf=>"Life",
                :cfi=>"FixCap",
                :ucap=>"unit_cap",
                :type=>"type",
                :xm=>"Max_capacity",
                :location=>"location")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    Pnames = convert(Array{String,1},sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"][:,1])
    locationcol = findfirst(col_names .== ndict[:location])[2]
    location = convert(Array{Int64,1},sheet[XLSX.encode_column_number(locationcol)*"3"*":"*XLSX.encode_column_number(locationcol)*"$stpr"][:,1])
    typecol = findfirst(col_names .== ndict[:type])[2]
    type = sheet[XLSX.encode_column_number(typecol)*"3"*":"*XLSX.encode_column_number(typecol)*"$stpr"][:,1]

    L=String[];G=String[];RS=String[];RW=String[];SE=String[];EB=String[];Dcov=String[];W=String[];C=String[];GB=String[];Dsep=String[];
    PO=String[];PWI=String[];PWL=String[];E=String[];F=String[];SHy=String[];B=String[];GI=String[];
    for i in 1:length(Pnames)
        if type[i]=="line"
            push!(L,Pnames[i])
        elseif type[i]=="OCGT"
            push!(G,Pnames[i])
        elseif type[i]=="solar"
            push!(RS,Pnames[i])
        elseif type[i]=="wind"
            push!(RW,Pnames[i])
        elseif type[i]=="battery"
            push!(SE,Pnames[i])
        elseif type[i]=="E_boiler"
            push!(EB,Pnames[i])
        elseif type[i]=="converter"
            push!(Dcov,Pnames[i])
        elseif type[i]=="NG_well"
            push!(W,Pnames[i])
        elseif type[i]=="compressor"
            push!(C,Pnames[i])
        elseif type[i]=="NG_boiler"
            push!(GB,Pnames[i])
        elseif type[i]=="separator"
            push!(Dsep,Pnames[i])
        elseif type[i]=="oil_exp_pump"
            push!(PO,Pnames[i])
        elseif type[i]=="water_inj_pump"
            push!(PWI,Pnames[i])
        elseif type[i]=="water_lift_pump"
            push!(PWL,Pnames[i])
        elseif type[i]=="gas_inj_compressor"
            push!(GI,Pnames[i])
        elseif type[i]=="electrolyser"
            push!(E,Pnames[i])
        elseif type[i]=="fuel_cell"
            push!(F,Pnames[i])
        elseif type[i]=="hydrogen_store"
            push!(SHy,Pnames[i])
        elseif type[i]=="bus"
            push!(B,Pnames[i])
    end
    end
    setfield!(ps,:L,L);setfield!(ps,:G,G);setfield!(ps,:RS,RS);setfield!(ps,:RW,RW);setfield!(ps,:SE,SE);setfield!(ps,:EB,EB);
    setfield!(ps,:Dcov,Dcov);setfield!(ps,:W,W);setfield!(ps,:C,C);setfield!(ps,:GB,GB);setfield!(ps,:Dsep,Dsep);
    setfield!(ps,:PO,PO);setfield!(ps,:PWI,PWI);setfield!(ps,:PWL,PWL);setfield!(ps,:E,E);setfield!(ps,:F,F);
    setfield!(ps,:SHy,SHy);setfield!(ps,:B,B);setfield!(ps,:GI,GI);

    PP=String[];PH=String[];PO=String[];L=String[];P=String[];
    for l in 1:length(location)
        if location[l]==0
            push!(PP,Pnames[l])
            push!(P,Pnames[l])
        elseif location[l]==1
            push!(PH,Pnames[l])
            push!(P,Pnames[l])
        elseif location[l]==2
            push!(PO,Pnames[l])
            push!(P,Pnames[l])
        elseif location[l]==3
            push!(L,Pnames[l])
            push!(P,Pnames[l])
        end
    end
    setfield!(ms,:PP,PP)
    setfield!(ms,:PH,PH)
    setfield!(ms,:PO,PO)
    setfield!(ms,:L,L)
    setfield!(ms,:P,P)
    setfield!(ps,:P,P)
    map_ppp=Dict{String,Int64}();map_php=Dict{String,Int64}();map_pop=Dict{String,Int64}();map_ll=Dict{String,Int64}()
    for pp in 1:length(PP)
        map_ppp[PP[pp]]=pp
    end
    for ph in 1:length(PH)
        map_php[PH[ph]]=length(PP)+ph
    end
    for po in 1:length(PO)
        map_pop[PO[po]]=length(PP)+length(PH)+po
    end
    for l in 1:length(L)
        map_ll[L[l]]=length(PP)+length(PH)+length(PO)+l
    end
    setfield!(mp,:map_ppp,map_ppp)
    setfield!(mp,:map_php,map_php)
    setfield!(mp,:map_pop,map_pop)
    setfield!(mp,:map_ll,map_ll)
    cicol = findfirst(col_names .== ndict[:ci])[2]
    cival = convert(Array{Float64,1},sheet[XLSX.encode_column_number(cicol)*"3"*":"*XLSX.encode_column_number(cicol)*"$stpr"][:,1])
    cfcol = findfirst(col_names .== ndict[:cf])[2]
    cfval = convert(Array{Float64,1},sheet[XLSX.encode_column_number(cfcol)*"3"*":"*XLSX.encode_column_number(cfcol)*"$stpr"][:,1])
    xmcol = findfirst(col_names .== ndict[:xm])[2]
    xmval = convert(Array{Float64,1},sheet[XLSX.encode_column_number(xmcol)*"3"*":"*XLSX.encode_column_number(xmcol)*"$stpr"][:,1])
    lfcol = findfirst(col_names .== ndict[:lf])[2]
    lfval = convert(Array{Float64,1},sheet[XLSX.encode_column_number(lfcol)*"3"*":"*XLSX.encode_column_number(lfcol)*"$stpr"][:,1])
    c_inv_0 = cival./lfval
    ucapcol = findfirst(col_names .== ndict[:ucap])[2]
    ucapval = convert(Array{Float64,1},sheet[XLSX.encode_column_number(ucapcol)*"3"*":"*XLSX.encode_column_number(ucapcol)*"$stpr"][:,1])
    cficol = findfirst(col_names .== ndict[:cfi])[2]
    cfival = convert(Array{Float64,1},sheet[XLSX.encode_column_number(cficol)*"3"*":"*XLSX.encode_column_number(cficol)*"$stpr"][:,1])


    ci = hcat(repeat(c_inv_0,1,NI1)).*10^3
    cf = cfval.*10^3
    xm0 = xmval.*10^-3
    setfield!(mp,:ci,ci)
    setfield!(mp,:cf,cf)
    setfield!(mp,:xm0,xm0)
    cfi = cfival./lfval
    xucap = ucapval.*10^-3
    setfield!(mp,:cfi,cfi)
    setfield!(mp,:xucap,xucap)

    udict = Dict(:zn =>"node",
                :type=>"node_type",
                :x0 =>"capacity_at_0",
                :x5 =>"capacity_at_5")  #column name in each sheet

    map_pp=Dict{String,Int64}();map_ph=Dict{String,Int64}();map_po=Dict{String,Int64}();map_l=Dict{String,Int64}()
    for pp in 1:length(PP)
        map_pp[PP[pp]]=pp
    end
    for ph in 1:length(PH)
        map_ph[PH[ph]]=ph
    end
    for po in 1:length(PO)
        map_po[PO[po]]=po
    end
    for l in 1:length(L)
        map_l[L[l]]=l
    end
    setfield!(mp,:map_pp,map_pp)
    setfield!(mp,:map_ph,map_ph)
    setfield!(mp,:map_po,map_po)
    setfield!(mp,:map_l,map_l)
    setfield!(pp,:map_l,map_l)

    xhpp=zeros(length(ms.PP),length(ms.ZP),NI2)
    for p in PP
        sheet = XFI[p]
        stt,stp = sheet.dimension.start,sheet.dimension.stop
        sttr,sttc = stt.row_number,stt.column_number
        stpr,stpc = stp.row_number,stp.column_number
        col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
        x5_col = findfirst(col_names .== udict[:x5])[2]
        x5_val = convert(Array{Float64,1},sheet[XLSX.encode_column_number(x5_col)*"3"*":"*XLSX.encode_column_number(x5_col)*"$stpr"][:,1])
        xhpp[map_pp[p],:,:] = hcat(repeat(x5_val,1,NI2)).*10^-3
    end
    setfield!(mp,:xhpp,xhpp)

    xhph=zeros(length(ms.PH),length(ms.ZH),NI2)
    for p in PH
        sheet = XFI[p]
        stt,stp = sheet.dimension.start,sheet.dimension.stop
        sttr,sttc = stt.row_number,stt.column_number
        stpr,stpc = stp.row_number,stp.column_number
        col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
        x5_col = findfirst(col_names .== udict[:x5])[2]
        x5_val = convert(Array{Float64,1},sheet[XLSX.encode_column_number(x5_col)*"3"*":"*XLSX.encode_column_number(x5_col)*"$stpr"][:,1])
        xhph[map_ph[p],:,:] = hcat(repeat(x5_val,1,NI2)).*10^-3
    end
    setfield!(mp,:xhph,xhph)

    xhpo=zeros(length(ms.PO),length(ms.ZO),NI2)
    for p in PO
        sheet = XFI[p]
        stt,stp = sheet.dimension.start,sheet.dimension.stop
        sttr,sttc = stt.row_number,stt.column_number
        stpr,stpc = stp.row_number,stp.column_number
        col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
        x5_col = findfirst(col_names .== udict[:x5])[2]
        x5_val = convert(Array{Float64,1},sheet[XLSX.encode_column_number(x5_col)*"3"*":"*XLSX.encode_column_number(x5_col)*"$stpr"][:,1])
        xhpo[map_po[p],:,:] = hcat(repeat(x5_val,1,NI2)).*10^-3
    end
    setfield!(mp,:xhpo,xhpo)
    xhl=zeros(length(ms.L),NI2)

    ############ read lines #################
    sheet = XFI["lines"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    x5_col = findfirst(col_names .== udict[:x5])[2]
    x5_val = convert(Array{Float64,1},sheet[XLSX.encode_column_number(x5_col)*"3"*":"*XLSX.encode_column_number(x5_col)*"$stpr"][:,1])
    xhl = hcat(repeat(x5_val,1,NI2)).*10^-3
    setfield!(mp,:xhl,xhl)
    setfield!(mp,:κ,5.)

    tx=Vector{Tuple{String,String}}()
    for i in 1:length(PP)
        for j in 1:length(ZP)
            push!(tx,(PP[i],ZP[j]))
        end
    end
    for i in 1:length(PH)
        for j in 1:length(ZH)
            push!(tx,(PH[i],ZH[j]))
        end
    end
    for i in 1:length(PO)
        for j in 1:length(ZO)
            push!(tx,(PO[i],ZO[j]))
        end
    end
    for l in 1:length(L)
        push!(tx,(L[l],L[l]))
    end
    for i in 1:length(th)
        push!(tx,(th[i],"-"))
    end
    setfield!(unc,:tx,tx)

    ############## gas turbine ############
    sheet = XFO_Electricity["gas_turbine"]
    ndict = Dict(:CG=>"VarOM",
                :CFuel=>"FuelCost",
                :EG=>"CO2_emissions",
                :ηG=>"efficiency",
                :GGR=>"ramping")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    Gnames = sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(Float64,sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1])
        setfield!(pp,vname,vvalue)
    end

    ############ battery ############
    sheet = XFO_Electricity["storage"]

    ndict = Dict(:ηSE=>"efficiency",
                :ΗSE=>"power_ratio",
                :CSE=>"VarOM")

    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    Bnames = sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(Float64,sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1])
        setfield!(pp,vname,vvalue)
    end

    # ############ time series ############
    sheet = XFO_Electricity["time_series"]

    ndict = Dict(1=>"start",
                2=>"end")

    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    Snames = sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"]
    setfield!(ps,:S,1:length(Snames))
    Hs = Array{Int64,2}(undef,maximum(ps.S),2)
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        Hs[:,vname] = convert(Array{Int64,1},sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][:,1])
    end
    setfield!(ps,:Hs,[Hs[s,1]:Hs[s,2] for s in ps.S])
    setfield!(ps,:H,1:8760)
    sH=Int64[];
    for s in ps.S
        for i in ps.Hs[s]
            push!(sH,i)
        end
    end
    setfield!(ps,:sH,sH)


    ############ renewable ############
    sheet = XFO_Electricity["wind"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    wvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])
    setfield!(pp,:RPW,wvalue)
    sheet = XFO_Electricity["solar"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    svalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])
    setfield!(pp,:RPS,svalue)
    ############ onshore electricity price ############
    sheet = XFO_Electricity["electricity_price"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    pricevalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])
    setfield!(pp,:CZO,pricevalue)
    ############ transmission line efficiency #########
    sheet = XFO_Electricity["line_loss"]
    ndict = Dict(:ηL=>"loss_factor")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    Gnames = sheet[XLSX.encode_column_number(sttc)*"3"*":"*XLSX.encode_column_number(sttc)*"$stpr"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(Array{Float64,1},sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][:,1])
        setfield!(pp,vname,vvalue)
    end
    ############ bus line incidence matrix ############
    sheet = XFO_Electricity["bus_line_matrix"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    A = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+1)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])
    setfield!(pp,:AE,A)

    ############ other electricity ############
    sheet = XFO_Electricity["other"]
    ndict = Dict(:CLShed => "load_shedding",
                :CGShed => "generation_shedding",
                :WS  =>"alpha_adjustment",
                :Ht=>"hour",
                :μD=>"demand_scaling",
                :μE=>"CO2_limit_scaling",
                :CCO2=>"CO2_tax",
                :ΜE=>"CO2_limit",
                :σRes=>"spinning_reserve")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ############ abatement cost ############
    sheet = XFO_Electricity["abatement_cost"]
    ndict = Dict(:CZOE => "abatement_cost_emission",
                :CZOP => "abatement_cost_power")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(Array{Float64,1},sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][:,1])
        setfield!(pp,vname,vvalue)
    end

    ############ natural gas ############
    sheet = XFO_NG["production"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*0.84
    setfield!(pp,:VD,Dvalue)
    sheet = XFO_NG["gas_injection"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*0.84
    setfield!(pp,:VGI,Dvalue)
    ############ gas well ############
    sheet = XFO_NG["well"]
    ndict = Dict(:CW =>"VarOM")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ############ compressor ############
    sheet = XFO_NG["compressor"]
    ndict = Dict(:γC =>"compression ratio",
                 :ηC => "efficiency",
                 :α => "polytropic exponent",
                 :CC => "VarOM")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ############ separator ############
    sheet = XFO_NG["separator"]
    ndict = Dict(:CSEP =>"VarOM",
                :κSEP=>"electricity_demand")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ############ NG_boiler ############
    sheet = XFO_NG["NG_boiler"]
    ndict = Dict(:CGB =>"VarOM",
                 :ηGB =>"efficiency")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ############ other ############
    sheet = XFO_NG["other"]
    ndict = Dict(:CNGShed =>"load_shed_cost",
                 :θNG => "gas_energy_content")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    ###### oil data ##########
    sheet = XFO_O["production"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*900
    setfield!(pp,:VOD,Dvalue)
    sheet = XFO_O["oil_pump"]
    ndict = Dict(:κPO=>"electricity_demand_po")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ###### water data ##########
    sheet = XFO_W["water_bore"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*1000
    setfield!(pp,:VWB,Dvalue)
    sheet = XFO_W["water_injection"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*1000
    setfield!(pp,:VWI,Dvalue)
    sheet = XFO_W["water_lift"]
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    Dvalue = convert(Array{Float64,2},sheet[XLSX.encode_column_number(sttc+1)*"$(sttr+2)"*":"*XLSX.encode_column_number(stpc)*"$stpr"])*1000
    setfield!(pp,:VWL,Dvalue)
    sheet = XFO_W["pump_wi"]
    ndict = Dict(:κPWI=>"electricity_demand_pwi")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end
    setfield!(pp,:VWL,Dvalue)
    sheet = XFO_W["pump_wl"]
    ndict = Dict(:κPWL=>"electricity_demand_pwl")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ############ heat data ############
    sheet = XFO_heat["efficiency"]
    ndict = Dict(:ηGh =>"heat_recovery_gas_turbine",
                :ηFh  =>"heat_recovery_fuel_cell",
                :ρSEP=>"heat_demand_separator",
                :ηEB=>"electricity_boiler_efficiency",
                :CHLShed=>"heat_shed_cost")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ############ electrolyser ############
    sheet = XFO_hydrogen["electrolyser"]
    ndict = Dict(:CE =>"VarOM",
                :ηES =>"efficiency_ES",
                :ηEF =>"efficiency_EF",
                :σE => "min_power")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ############ hydrogen storage ############
    sheet = XFO_hydrogen["hydrogen_store"]
    ndict = Dict(:CSHy =>"VarOM",
                :σSHy =>"cushion_gas")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ############ fuel cell ############
    sheet = XFO_hydrogen["fuel_cell"]
    ndict = Dict(:CF =>"VarOM",
                 :ηF => "efficiency",
                 :FR =>"ramping")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ############ other hydrogen data ############
    sheet = XFO_hydrogen["other"]
    ndict = Dict(:θHY=>"hydrogen_energy_content")
    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]
    for vname in keys(ndict)
        vcol = findfirst(col_names .== ndict[vname])[2]
        vvalue = convert(typeof(getfield(pp,vname)),sheet[XLSX.encode_column_number(vcol)*"3"*":"*XLSX.encode_column_number(vcol)*"$stpr"][1,1])
        setfield!(pp,vname,vvalue)
    end

    ########## settings #########
    sheet = XFA["settings"]

    ndict = Dict(:Imax=>"Max_Iteration",
                 :tol =>"Convergence_Tol",
                 :alg =>"Type",
                 :γs  =>"stabilization_factor",
                 :w   =>"Subproblems_each_iter")

    stt,stp = sheet.dimension.start,sheet.dimension.stop
    sttr,sttc = stt.row_number,stt.column_number
    stpr,stpc = stp.row_number,stp.column_number
    col_names = sheet[XLSX.encode_column_number(sttc)*"1"*":"*XLSX.encode_column_number(stpc)*"1"]

    Icol = findfirst(col_names .== ndict[:Imax])[2]
    J = convert(Int64,sheet[XLSX.encode_column_number(Icol)*"3"*":"*XLSX.encode_column_number(Icol)*"$stpr"][1,1])
    ϵcol = findfirst(col_names .== ndict[:tol])[2]
    ϵ = convert(Float64,sheet[XLSX.encode_column_number(ϵcol)*"3"*":"*XLSX.encode_column_number(ϵcol)*"$stpr"][1,1])
    acol = findfirst(col_names .== ndict[:alg])[2]
    alg = convert(Int64,sheet[XLSX.encode_column_number(acol)*"3"*":"*XLSX.encode_column_number(acol)*"$stpr"][1,1])
    γcol = findfirst(col_names .== ndict[:γs])[2]
    γs = convert(Float64,sheet[XLSX.encode_column_number(γcol)*"3"*":"*XLSX.encode_column_number(γcol)*"$stpr"][1,1])
    wcol = findfirst(col_names .== ndict[:w])[2]
    w = convert(Int64,sheet[XLSX.encode_column_number(wcol)*"3"*":"*XLSX.encode_column_number(wcol)*"$stpr"][1,1])

    return J,ϵ,alg,γs,unc,ms,mp,ps,pp
end
