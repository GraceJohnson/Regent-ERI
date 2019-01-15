import "regent"

-- Helper module to handle command line arguments
local ERIConfig = require("eri_config")

local c = regentlib.c
local cstring = terralib.includec("string.h")
local cmath = terralib.includec("math.h")

local sqrt = regentlib.sqrt(double)
local pow = regentlib.pow(double)
local erf = regentlib.erf(double)

MAX_BASIS_SIZE = 10
MAX_NUM_ATOMS = 1000
MAX_NUM_PRIMS = MAX_BASIS_SIZE * MAX_NUM_ATOMS

--------------------------
---- Field Space Defs ----
--------------------------

fspace Atom {
  x  :  double;
  y  :  double;
  z  :  double;
}

fspace Primitive {
  x  :  double;
  y  :  double;
  z  :  double;
  w  :  double; -- primitive coefficient
  c  :  double; -- normalized primitive coefficient
  e  :  double; -- exponent of Gaussian
  AOid : int1d; -- used to map/partition/color the region of primitives into atomic orbitals
}

fspace Pairlist {
  -- see Boys formula from 2008 TC paper for details on these attributes
  Eij       :  double; -- ei + ej
  Kij       :  double; -- normalized exponential 
  Rijx      :  double; -- normalized distance
  Rijy      :  double; -- normalized distance
  Rijz      :  double; -- normalized distance
  AOpair_id :  int1d; -- used to map/partition/color the region of primitives into atomic orbitals
  xpart_id  :  int1d;
  ypart_id  :  int1d;
}

-- like CSV storage
fspace AOpairlist {
  start_idx : int;
  end_idx   : int;
}

fspace SIntegral {
  s_int  :  double;
}

fspace EW_data {
  ew    :   double;
}

fspace AOidx_data {
  aoidx    :   int;
}

---------------------------
------ UTILITY TASKS ------
---------------------------

terra skip_header(f : &c.FILE)
  c.fscanf(f, "%*llu\n%*[^\n]\n")
end

terra skip_basis_header(f : &c.FILE)
  var line : int8[512]
  c.fscanf(f,"%[^\n]\n", &line)
  if (cstring.strcmp(line, "****") == 0) then return true
  else return false end
end

-- Returns false if there was an error reading, e.g. if something besides H
terra read_xyz(f : &c.FILE, xyz : &double)
  var atom_type : int8[512]
  c.fscanf(f, "%s%lf%lf%lf\n", &atom_type, &xyz[0], &xyz[1], &xyz[2])
  if not (cstring.strcmp(atom_type, "H") == 0) then return false 
  else return true end
end

terra read_basis_file(f       : &c.FILE,
                      EW      : &double,
                      AOidx   : &int,
                      counter : &int8)
  var mobj1 : int8[512]
  var mobj2 : int8[512]
  c.fscanf(f,"%s%s%*[^\n]\n", &mobj1, &mobj2)
  if (cstring.strcmp(mobj1, "****") == 0) or (cstring.strcmp(mobj1, "He") == 0) then  
    return true
  elseif (cstring.strcmp(&mobj1[0], "H") == 0) then
    return false 
  elseif (cstring.strcmp(&mobj1[0], "S") == 0) then
    var shell_done = false
    AOidx[counter[1]] = counter[0]
    counter[1] = counter[1] + 1
    while not shell_done do
      var mobj3 : double[1] 
      var mobj4 : double[1]
      mobj3[0] = 0
      mobj4[0] = 0
      var shell = c.fscanf(f, "%lf%lf\n", &mobj3, &mobj4)
      
      if not(shell == 2) then shell_done = true 
      else
        EW[counter[0]] = double(mobj3[0])
        EW[counter[0]+1] = double(mobj4[0])
        counter[0] = counter[0] + 2
      end
    end
    return false 
  else
    c.printf("ERROR: Error reading basis file")
    return true
  end
end


task factorize(parallelism : int) : int2d
  var limit = [int](cmath.sqrt([double](parallelism)))
  var size_x = 1
  var size_y = parallelism
  for i = 1, limit + 1 do
    if parallelism % i == 0 then
      size_x, size_y = i, parallelism / i
      if size_x > size_y then
        size_x, size_y = size_y, size_x
      end
    end
  end
  return int2d { size_x, size_y }
end


task write_output(r_eri     :  region(ispace(int2d), SIntegral),
                  filename  :  int8[512], 
                  M         :  int)
where 
  reads(r_eri.s_int)
do
  var f = c.fopen(filename, "w")
  for i = 0, M do
    for j = 0, M do
      c.fprintf(f, "%d  %d  %.10f\n", i, j, r_eri[{i, j}].s_int)
    end
  end
  c.fclose(f)
  

end


----------------------------
---- ERI IMPLEMENTATION ----
----------------------------

task initialize_system(r_atoms     : region(ispace(int1d), Atom),
                       r_EW        : region(ispace(int1d), EW_data),
                       r_AOidx     : region(ispace(int1d), AOidx_data),
                       geom_file   : int8[512],
                       basis_file  : int8[512])

where 
  reads writes(r_atoms, r_EW, r_AOidx)
do
  var ts_start = c.legion_get_current_time_in_micros()

  var EW : double[MAX_BASIS_SIZE * 2] -- stores raw data from basis file
  var AOidx: int[MAX_BASIS_SIZE] -- stores indicies where each new contraction starts in EW

  -- Read in input geometry from XYZ file
  var f = c.fopen(geom_file, "rb")
  skip_header(f)
  var xyz : double[3]
  for atom in r_atoms do
    --regentlib.assert(read_xyz(f, xyz), "ERI implemented only for H atoms! \
    --              Please provide geom with only H atoms.\n")
    read_xyz(f,xyz)
    atom.x = xyz[0]
    atom.y = xyz[1]
    atom.z = xyz[2]
  end
  c.fclose(f)

  -- Read in basis set info (for H atoms only)
  var fb = c.fopen(basis_file, "rb")
  var done = false
  while not done do
    done = skip_basis_header(fb)
  end
  done = false
  var counter : int8[2]
  counter[0] = 0 -- the coutner for EW array
  counter[1] = 0 -- the counter for AOidx array
  while not done do
    done = read_basis_file(fb, EW, AOidx, counter)
  end
  c.fclose(fb)

  -- copy array values to regions
  for i in r_EW.ispace do
    r_EW[i].ew = EW[int(i)]
  end
  for i in r_AOidx.ispace do
    r_AOidx[i].aoidx = AOidx[int(i)]
  end
  

  var ts_stop = c.legion_get_current_time_in_micros()
  c.printf("System initialization took %.4f sec\n", (ts_stop - ts_start) * 1e-6)

  return counter
end


task initialize_prims(r_atoms      : region(ispace(int1d), Atom),
                      r_prims      : region(ispace(int1d), Primitive),
                      r_EW        : region(ispace(int1d), EW_data),
                      r_AOidx     : region(ispace(int1d), AOidx_data),
                      num_prims    : int,
                      num_AO       : int,
                      num_atoms    : uint64)
where 
  reads (r_atoms, r_EW, r_AOidx),
  reads writes (r_prims) 
do
  var ts_start = c.legion_get_current_time_in_micros()

  -- compute normalized coefficients for contracted orbitals
  ----------------------------------------------------------------------------
  var C : double[MAX_NUM_PRIMS]
  var count = 0
  -- NOTE: this is only for L = 0 case (angular momentum = 0 for s orbitals)
  var pi32 = cmath.pow(cmath.M_PI, 1.5)

  -- loop over all atomic orbitals
  for i = 0, num_AO do
    -- loop over primitives within AO to normalize prim coeffs
    for j = r_AOidx[i].aoidx/2, r_AOidx[i+1].aoidx/2 do 
      C[j] = r_EW[2*j+1].ew * c.sqrt(c.pow(2.0 * r_EW[2*j].ew, 1.5) / pi32)
    end
    -- normalize for overlap of one basis function with another (pair) -- unnecessary for s?
    var V = 0.0
    for j1 = r_AOidx[i].aoidx/2, r_AOidx[i+1].aoidx/2 do 
      for j2 = r_AOidx[i].aoidx/2, r_AOidx[i+1].aoidx/2 do 
        V += C[j1] * C[j2] * cmath.pow(r_EW[2*j1].ew + r_EW[2*j2].ew, -1.5) * pi32
      end
    end
    V = cmath.pow(V, -1.0/2.0) -- V just ends up being 1.0
    for j = r_AOidx[i].aoidx/2, r_AOidx[i+1].aoidx/2 do 
      C[j] *= V
    end
  end

  ----------------------------------------------------------------------------

  var basis_count = 0
  var atom_count = 0
  
  -- loop through atoms and populate prims array
  var geoms : double[MAX_NUM_PRIMS * 3]
  var prim_count = 0
  for atom in r_atoms do
    for i = 0, int(num_prims) do 
      geoms[prim_count] = atom.x 
      geoms[prim_count + 1] = atom.y 
      geoms[prim_count + 2] = atom.z 
      prim_count = prim_count + 3
    end
  end

  -- populate prim region from prim array
  var geom_idx = 0
  var basis_idx = 0
  for prim in r_prims do
    prim.x = geoms[geom_idx]
    prim.y = geoms[geom_idx + 1]
    prim.z = geoms[geom_idx + 2]
    prim.e = r_EW[basis_idx].ew
    prim.w = r_EW[basis_idx + 1].ew
    prim.c = C[basis_idx/2]

    geom_idx = geom_idx + 3
    basis_idx = basis_idx + 2
    if basis_idx >= num_prims*2 then
      basis_idx = 0
    end
  end

  -- add mapping from primitive to AO basis
  var ao_id = 0
  var ao_map : int[MAX_NUM_PRIMS]
  var idx = 0
  while ao_id < num_AO * num_atoms do
    for i = 0, num_AO do
      for j = r_AOidx[i].aoidx/2, r_AOidx[i+1].aoidx/2 do 
        ao_map[idx] = ao_id
        idx += 1
      end
      ao_id += 1
    end
  end
  
  for i in r_prims.ispace do
    r_prims[i].AOid = ao_map[int(i)]
  end

  var ts_stop = c.legion_get_current_time_in_micros()
  c.printf("Initialization of primitives took %.4f sec\n", (ts_stop - ts_start) * 1e-6)

end


task compute_pairlists(r_prims      : region(ispace(int1d), Primitive),
                       r_pairlists  : region(ispace(int1d), Pairlist),
                       N            : int)
where 
  reads (r_prims), 
  reads writes (r_pairlists)
do
  var ts_start = c.legion_get_current_time_in_micros()

  -- loop over AO, and make pairs within each AO
  -- just add to pairlist region, include AOid, then partition later 
  var p_AO = partition(r_prims.AOid, ispace(int1d, N))

  var count = 0
  var pair_id = 0
  for AO1 in p_AO.colors do
    for AO2 in p_AO.colors do
      if AO1 <= AO2 then
        for i1 in p_AO[AO1].ispace do
          for i2 in p_AO[AO2].ispace do
            -- Fill out attributes of Pairlist
            var e1 = r_prims[i1].e
            var e2 = r_prims[i2].e
            var c1 = r_prims[i1].c
            var c2 = r_prims[i2].c
            var x1 = r_prims[i1].x
            var y1 = r_prims[i1].y
            var z1 = r_prims[i1].z
            var x2 = r_prims[i2].x
            var y2 = r_prims[i2].y
            var z2 = r_prims[i2].z
            -- Euclidian distance squared
            var dist = cmath.pow(x2-x1, 2.0) + cmath.pow(y2-y1, 2.0) + cmath.pow(z2-z1, 2.0) 
            var exp = cmath.exp(-e1*e2/(e1 + e2) * dist)
            var einv = 1.0/(e1 + e2) 
            r_pairlists[count].Eij = e1 + e2
            --r_pairlists[count].Kij = c1*c2 * einv * exp   -- this was part of the error  
            r_pairlists[count].Kij = c1 * c2 * exp  
            r_pairlists[count].Rijx = einv * (e1*x1 + e2*x2)  
            r_pairlists[count].Rijy = einv * (e1*y1 + e2*y2)
            r_pairlists[count].Rijz = einv * (e1*z1 + e2*z2)
            r_pairlists[count].AOpair_id = pair_id
            count += 1
          end
        end
        pair_id += 1
      end
    end
  end
  
  var ts_stop = c.legion_get_current_time_in_micros()
  c.printf("Computing pairlists took %.4f sec\n", (ts_stop - ts_start) * 1e-6)

end


__demand(__cuda)
task compute_integrals(r_eri          : region(ispace(int2d), SIntegral),
                       rx_pairlists   : region(ispace(int1d), Pairlist),
                       ry_pairlists   : region(ispace(int1d), Pairlist),
                       rx_AOpairlists : region(ispace(int1d), AOpairlist),
                       ry_AOpairlists : region(ispace(int1d), AOpairlist))
where 
  reads (rx_pairlists, ry_pairlists, rx_AOpairlists, ry_AOpairlists),
  writes (r_eri.s_int)
do
  
  var pi12 = c.sqrt(c.M_PI)
  var pi52 = c.pow(c.M_PI, 5.0/2.0)

  var xsize = rx_AOpairlists.ispace.volume 
  var ysize = ry_AOpairlists.ispace.volume 

  var AO2d_ispace = ispace(int1d, xsize*ysize)
  
  for AO2d in AO2d_ispace do
    var x = (AO2d % xsize) + rx_AOpairlists.ispace.bounds.lo
    var y = (AO2d / xsize) + ry_AOpairlists.ispace.bounds.lo
    var AO1 = rx_AOpairlists[x]
    var AO2 = ry_AOpairlists[y]

    var integral = 0.0
    for i1 = AO1.start_idx, AO1.end_idx do
      var x1 = rx_pairlists[i1].Rijx
      var y1 = rx_pairlists[i1].Rijy
      var z1 = rx_pairlists[i1].Rijz
      var E1 = rx_pairlists[i1].Eij
      var K1 = rx_pairlists[i1].Kij
      for i2 = AO2.start_idx, AO2.end_idx do
        var x2 = ry_pairlists[i2].Rijx
        var y2 = ry_pairlists[i2].Rijy
        var z2 = ry_pairlists[i2].Rijz
        var E2 = ry_pairlists[i2].Eij
        var K2 = ry_pairlists[i2].Kij
        var dist = pow(x2 - x1, 2.0) + pow(y2 - y1, 2.0) + pow(z2 - z1, 2.0)
        var t = sqrt(E1 * E2/(E1 + E2) * dist)
        var Einv = 1.0/(E1 * E2 * sqrt(E1 + E2))
        if t < 1.0e-10 then -- erf(sqrt(0))/sqrt(0) = 1
          integral += 2.0 * pi52 * Einv * K1 * K2 
        else
          integral += 2.0 * pi52 * Einv * K1 * K2 * (pi12/2.0) * (1.0/t) * erf(t) 
        end
      end
    end 
    r_eri[{x, y}].s_int = integral
  end

end


task check_normalization(r_prims  : region(ispace(int1d), Primitive),
                         N        : int)
where 
  reads (r_prims) 
do
  var p_AO = partition(r_prims.AOid, ispace(int1d, N))

  c.printf("Printing overlaps...\n\n")
  var count = 0
  for AO in p_AO.colors do
    var overlap = 0.0
    for i1 in p_AO[AO].ispace do
      for i2 in p_AO[AO].ispace do
        var e1 = r_prims[i1].e
        var e2 = r_prims[i2].e
        var c1 = r_prims[i1].c
        var c2 = r_prims[i2].c
        var oe = 1.0/(e1 + e2)
        overlap += c1 * c2 * c.pow(cmath.M_PI * oe, 1.5)
      end
    end
    c.printf("AO : %d    overlap : %.6f\n", AO, overlap)
  end

end





------------------------
---- TOP LEVEL TASK ----
------------------------
 
task toplevel()
  var config : ERIConfig
  config:initialize_from_command()
  c.printf("*********************************************\n")
  c.printf("*      Two-Electron Repulsion Integrals     *\n")
  c.printf("*                                           *\n")
  c.printf("* Input Geom       : %20s   *\n",  config.geom)
  c.printf("* Number of Atoms  : %20lu   *\n",  config.num_atoms)
  c.printf("* Basis Set        : %20s   *\n",  config.basis_name)
  c.printf("* # Parallel Tasks : %20u   *\n",   config.parallelism)
  c.printf("*********************************************\n")

  var tot_ts_start = c.legion_get_current_time_in_micros()

  -- Create a region of atoms
  var r_atoms = region(ispace(int1d, config.num_atoms), Atom)

  -- Set up system (read in geom and basis)

  var r_EW = region(ispace(int1d, MAX_BASIS_SIZE*2), EW_data)
  var r_AOidx = region(ispace(int1d, MAX_BASIS_SIZE), AOidx_data)

  var counter: int8[2]
  counter = initialize_system(r_atoms, r_EW, r_AOidx, config.geom, config.basis)
  var EW_size = counter[0]
  var num_AO_per_atom = counter[1] -- number of atomic orbitals
  var N = config.num_atoms*num_AO_per_atom
  var M = N * (N + 1)/2
  r_AOidx[num_AO_per_atom].aoidx = EW_size

  -- Create a region of primitives of size num_prims
  var num_prims = config.num_atoms*(EW_size/2)
  var r_prims = region(ispace(int1d, num_prims), Primitive)

  -- load data into r_prims and normalize gaussian coefficients 
  initialize_prims(r_atoms, r_prims, r_EW, r_AOidx, EW_size/2, num_AO_per_atom, config.num_atoms)

  -- Create a region of pairlists
  var num_pairs = 0 -- count up total number of primitive pairs
  var prims : int[MAX_BASIS_SIZE]
  for i = 0, num_AO_per_atom do
    var n = r_AOidx[i+1].aoidx/2 - r_AOidx[i].aoidx/2
    prims[i] = n
  end
  for i = 0, N do
    for j = 0, N do
      if i <= j then
        num_pairs += prims[i % num_AO_per_atom] * prims[j % num_AO_per_atom]
      end
    end
  end
  var r_pairlists = region(ispace(int1d, num_pairs), Pairlist) 

  -- Compute pairlists (simple compared to 2ERI, no need to parallelize)
  compute_pairlists(r_prims, r_pairlists, N)

  -- create CSR-like array in AOpairlist to index into r_pairlists for atomic orbitals
  var r_AOpairlists = region(ispace(int1d, M), AOpairlist) 
  var curr_AO : int1d = 0
  var count = 0
  r_AOpairlists[0].start_idx = curr_AO
  for i in r_pairlists.ispace do
    if r_pairlists[i].AOpair_id == curr_AO + 1 then
      r_AOpairlists[count].end_idx = i
      count += 1
      if count < M then
        r_AOpairlists[count].start_idx = i
        curr_AO += 1
      end
    end
  end
  r_AOpairlists[M-1].end_idx = num_pairs

  -- Create a region to store computed two-electron repulsion integral
  var r_eri = region(ispace(int2d, {M, M}), SIntegral)

  -- block r_eri and r_pairlists according to config.parallelism
  var block2d = factorize(config.parallelism)
  
  -- partition r_eri and r_AOpairlist
  var p_eri = partition(equal, r_eri, ispace(int2d, block2d))
  var px_AOpairlists = partition(equal, r_AOpairlists, ispace(int1d, block2d.x)) 
  var py_AOpairlists = partition(equal, r_AOpairlists, ispace(int1d, block2d.y)) 

  -- label r_pairlist for parallelizing
  for color in px_AOpairlists.colors do
    var chunk = px_AOpairlists[color]
    for pair in chunk do
      for i = pair.start_idx, pair.end_idx do
        r_pairlists[i].xpart_id = color
      end
    end 
  end
  for color in py_AOpairlists.colors do
    var chunk = py_AOpairlists[color]
    for pair in chunk do
      for i = pair.start_idx, pair.end_idx do
        r_pairlists[i].ypart_id = color
      end
    end 
  end

  -- now we can partition the actual pairlist region based on label
  var px_pairlists = partition(r_pairlists.xpart_id, ispace(int1d, block2d.x))
  var py_pairlists = partition(r_pairlists.ypart_id, ispace(int1d, block2d.y))
  


  -- ERI Calculation ---------------

  __fence(__execution, __block) -- block for timing
  var ts_start = c.legion_get_current_time_in_micros()

  -- Loop over pairlists to compute integral
  for color in p_eri.colors do
    compute_integrals(p_eri[color], 
                      px_pairlists[color.x], 
                      py_pairlists[color.y],   
                      px_AOpairlists[color.x],
                      py_AOpairlists[color.y])
  end

  __fence(__execution, __block) -- block for timing
  var ts_stop = c.legion_get_current_time_in_micros()
  c.printf("ERI Calculation took %.4f sec\n", (ts_stop - ts_start)* 1e-6)


  var tot_ts_stop = c.legion_get_current_time_in_micros()
  c.printf("Total Calculation took %.4f sec\n", (tot_ts_stop - tot_ts_start)* 1e-6)

  -- write output ERI matrix to file
  if config.dump_output then write_output(r_eri, config.output, M) end

end

regentlib.start(toplevel)
