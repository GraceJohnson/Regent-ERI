import "regent"

local c = regentlib.c

local util = {}

struct ERIConfig
{
  geom          : int8[512],
  basis         : int8[512],
  output        : int8[512],
  dump_output   : bool,
  basis_name    : int8[512],
  num_atoms     : uint64;
  parallelism   : uint32;
}

local cstring = terralib.includec("string.h")

terra print_usage_and_abort()
  c.printf("Usage: regent eri.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h            : Print the usage and exit.\n")
  c.printf("  -i {file.xyz} : Use {file.xyz} as input geometry.\n")
  c.printf("  -b {file.gbs} : Use {file.gbs} as input basis set.\n")
  c.printf("  -o {file}     : Save output to {file}.\n")
  c.printf("  -p {value}    : Set the number of parallel tasks to {value}.\n")
  c.exit(0)
end

terra file_exists(filename: rawstring)
  var file = c.fopen(filename, "rb")
  if file == nil then return false end
  c.fclose(file)
  return true
end

terra ERIConfig:initialize_from_command()
  var geom_given = false
  var basis_given = false

  self.dump_output = false
  self.parallelism = 1

  var args = c.legion_runtime_get_input_args()
  var i = 1 
  while i < args.argc do
    if cstring.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()
    elseif cstring.strcmp(args.argv[i], "-i") == 0 then
      i = i + 1

      var file = c.fopen(args.argv[i], "rb")
      if file == nil then
        c.printf("File '%s' doesn't exist!\n", args.argv[i])
        c.abort()
      end
      cstring.strcpy(self.geom, args.argv[i])
      c.fscanf(file, "%llu\n", &self.num_atoms)
      geom_given = true
      c.fclose(file)
    elseif cstring.strcmp(args.argv[i], "-b") == 0 then
      i = i + 1

      var file = c.fopen(args.argv[i], "rb")
      if file == nil then
        c.printf("File '%s' doesn't exist!\n", args.argv[i])
        c.abort()
      end
      cstring.strcpy(self.basis, args.argv[i])
      -- capture basis name from filename (strip off suffix)
      c.fscanf(file, "%*s\n\n%*s\n%*s%s", &self.basis_name)
      --cstring.strncpy(self.basis_name, self.basis, (cstring.strlen(self.basis) - 4) )
      basis_given = true
      c.fclose(file)
    elseif cstring.strcmp(args.argv[i], "-o") == 0 then
      i = i + 1
      if file_exists(args.argv[i]) then
        c.printf("File '%s' already exists!\n", args.argv[i])
        c.abort()
      end
      cstring.strcpy(self.output, args.argv[i])
      self.dump_output = true
    elseif cstring.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.parallelism = c.atoi(args.argv[i])
    end
    i = i + 1
  end
  if not geom_given then
    c.printf("Input XYZ file must be given!\n\n")
    print_usage_and_abort()
  end
  if not basis_given then
    c.printf("Basis set file must be given!\n\n")
    print_usage_and_abort()
  end
end

return ERIConfig


