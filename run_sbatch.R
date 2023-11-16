# run_sbatch.R
cat('.................................................................\n')
cat('888d888.888..888.88888b..........................................\n')
cat('888P"...888..888.888."88b........................................\n')
cat('888.....888..888.888..888..................flowplyr v1.0.........\n')
cat('888.....Y88b.888.888..888........................................\n')
cat('888......"Y88888.888..888........................................\n')
cat('.................................................................\n')
cat('.................................................................\n')
cat('.................................................................\n')
cat('.........888...............888............888..........8888888b..\n')
cat('.........888...............888............888..........888...Y88b\n')
cat('.........888...............888............888..........888....888\n')
cat('.d8888b..88888b....8888b...888888..d8888b.88888b.......888...d88P\n')
cat('88K......888."88b....."88b.888...d88P"....888."88b.....8888888P".\n')
cat('"Y8888b..888..888..d888888.888...888......888..888.....888.T88b..\n')
cat('.....X88.888.d88P.888..888.Y88b..Y88b.....888..888.d8b.888..T88b.\n')
cat('..88888P.88888P".."Y888888.."Y888."Y8888P.888..888.Y8P.888...T88b\n')

parser <- argparser::arg_parser(
  description='run_sbatch.R - launch multiple shell scripts using sbatch')
parser <- argparser::add_argument(parser, arg="--folder", 
                                  type="character",
                                  help = "folder where bash shell scripts are")
parser <- argparser::add_argument(parser, arg="--pattern", 
                                  type="character",
                                  default = 'extract_.*.sh$',
                                  help = "pattern that script must match")
parser <- argparser::add_argument(parser, arg="--dry_run", 
                                  type="logical",
                                  default = TRUE,
                                  help = "Preview rather than actually kick of jobs")
params   <- argparser::parse_args(parser)
sbatch_shell_scripts = list.files(params$folder, pattern = params$pattern)
for (script in sbatch_shell_scripts){
  print( file.path(params$folder, script))
  if(params$dry_run){
    system(paste("more", file.path(params$folder, script)))
  }else{
    system(paste("sbatch", file.path(params$folder, script)))
  }
}
  
