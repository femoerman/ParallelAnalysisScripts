######################################################################
# R script for analysing video files with BEMOVI (www.bemovi.info)
#
# Felix Moerman
#
# July 2018
######################################################################
rm(list=ls())
setwd("/media/mendel-himself/LuisJ/Luis_stability_2018")

# load package
#library(devtools)
#install_github("efronhofer/bemovi", ref="experimental")
library(bemovi)
library(parallel)
library(doParallel)
library(foreach)

#Define memory to be allocated
mem.per.identifier <- c(30000)
mem.per.linker <- c(10000)
mem.per.overlay <- c(60000)
machine.ram <- c(250000)

# UNIX
# set paths to ImageJ and particle linker standalone
IJ.path.main <- "/home/mendel-himself/ImageJ"
to.particlelinker.main <- "/home/mendel-himself/ParticleLinker"

# directories and file names
video.description.folder <- "0_video_description/"
video.description.file <- "video_description.txt"
raw.video.folder <- "1_raw/"
particle.data.folder <- "2_particle_data/"
trajectory.data.folder <- "3_trajectory_data/"
temp.overlay.folder <- "4a_temp_overlays/"
overlay.folder <- "4_overlays/"
merged.data.folder <- "5_merged_data/"
ijmacs.folder <- "ijmacs/"

# video frame rate (in frames per second)
fps <- 25
# length of video (in frames)
total_frames <- 500

# measured volume (in microliter)
measured_volume <- 34.4 # for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
#measured_volume <- 14.9 # for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III

# size of a pixel (in micrometer)
pixel_to_scale <- 4.05 # for Leica M205 C with 1.6 fold magnification, sample height 0.5 mm and Hamamatsu Orca Flash 4
#pixel_to_scale <- 3.79 # for Nikon SMZ1500 with 2 fold magnification, sample height 0.5 mm and Canon 5D Mark III

# specify video file format (one of "avi","cxd","mov","tiff")
# bemovi only works with avi and cxd. other formats are reformated to avi below
video.format <- "cxd"

# setup
difference.lag <- 10
thresholds <- c(10,255) # don't change the second value
#thresholds <- c(50,255)

# MORE PARAMETERS (USUALLY NOT CHANGED)
######################################################################
# FILTERING PARAMETERS 
# optimized for Perfex Pro 10 stereomicrocope with Perfex SC38800 (IDS UI-3880LE-M-GL) camera
# tested stereomicroscopes: Perfex Pro 10, Nikon SMZ1500, Leica M205 C
# tested cameras: Perfex SC38800, Canon 5D Mark III, Hamamatsu Orca Flash 4
# tested species: Tet, Col, Pau, Pca, Eug, Chi, Ble, Ceph, Lox, Spi

# min and max size: area in pixels
particle_min_size <- 5
particle_max_size <- 1000

# number of adjacent frames to be considered for linking particles
trajectory_link_range <- 3
# maximum distance a particle can move between two frames
trajectory_displacement <- 16

# these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)
filter_min_net_disp <- 25
filter_min_duration <- 1
filter_detection_freq <- 0.1
filter_median_step_length <- 3

#From here on: difference with other script
#Get a list with all subfolders to be analysed
folderlist <- list.dirs(path = getwd(), full.names = TRUE, recursive = FALSE)

# setup
difference.lag <- 10
thresholds <- c(10,255) # don't change the second value

to.particlelinker <- to.particlelinker.main

######################################################################
#copy IJ to new folder, as well as identification script
#Modify particle identification script
ident = readLines("Parallelised_analysis_script_FM_27July2018_Identification.R")
ident_original <- readLines("Parallelised_analysis_script_FM_27July2018.R")

#Prepare code that needs to be changed in identification script
{
p1 <- "# video frame rate (in frames per second)"                                                                              
p2 <- paste("fps <-", fps)                                                                                                           
p3 <- "# length of video (in frames)"                                                                                          
p4 <- paste("total_frames <-", total_frames)                                                                                                    
p5 <- ""                                                                                                                    
p6 <- "# measured volume (in microliter)"                                                                                      
p7 <- ""
p8 <- paste("measured_volume <-",  measured_volume) 
p9 <- ""                                                                                                                       
p10 <- "# size of a pixel (in micrometer)"                                                                                      
p11 <- paste("pixel_to_scale <-", pixel_to_scale )
p12 <- ""      
p13 <- ""                                                                                                                       
p14 <- "# specify video file format (one of \"avi\",\"cxd\",\"mov\",\"tiff\")"                                                  
p15 <- "# bemovi only works with avi and cxd. other formats are reformated to avi below"                                        
p16 <- paste("video.format <-", video.format)                                                                                               
p17 <- ""                                                                                                                       
p18 <- "# setup"                                                                                                                
p19 <- paste("difference.lag <-", difference.lag)                                                                                                
p20 <- paste("thresholds <- c(" , thresholds[1],",", thresholds[2], ") # don't change the second value")                                                                
p21 <- ""
text1 <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21)

p1 <- "# min and max size: area in pixels"                                                                                                           
p2 <- paste("particle_min_size <-", particle_min_size)                                                                                                                      
p3 <- paste("particle_max_size <-", particle_max_size)                                                                                                                 
p4 <- ""                                                                                                                                             
p5 <- "# number of adjacent frames to be considered for linking particles"                                                                           
p6 <- paste("trajectory_link_range <-", trajectory_link_range)                                                                                                                   
p7 <- "# maximum distance a particle can move between two frames"                                                                                    
p8 <- paste("trajectory_displacement <-", trajectory_displacement)                                                                                                                
p9 <- ""                                                                                                                                             
p10 <- "# these values are in the units defined by the parameters above: fps (seconds), measured_volume (microliters) and pixel_to_scale (micometers)"
p11 <-  paste("filter_min_net_disp <-", filter_min_net_disp)                                                                                                                  
p12 <- paste("filter_min_duration <-", filter_min_duration)                                                                                                                      
p13 <- paste("filter_detection_freq <-", filter_detection_freq)                                                                                                                     
p14 <- paste("filter_median_step_length <-", filter_median_step_length)                                                                                                               
p15 <- ""                                                                                                                                             
p16 <- "#From here on: difference with other script"                                                                                                  
p17 <- "#Get a list with all subfolders to be analysed"                                                                                               
p18 <- "folderlist <- list.dirs(path = getwd(), full.names = TRUE, recursive = FALSE)"                                                                
p19 <- ""                                                                                                                                             
p20 <- "# setup"                                                                                                                                      
p21 <- paste("difference.lag <-", difference.lag)                                                                                                                         
p22 <- paste("thresholds <- c(" , thresholds[1],",", thresholds[2], ") # don't change the second value")  

text2 <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22)
}
ident[32:52] <- text1
ident[61:82] <- text2
ident[85] <- paste("mem.per.identifier <- c(", mem.per.identifier, ")")


for (foldername in folderlist){
  ident[87] <- paste("to.data <- '", foldername, "/'", sep="")
  writeLines(ident, paste(foldername, "/Parallelised_analysis_script_FM_27July2018_Identification.R", sep=""))
  newIJ <- paste(foldername, "/IJ", sep="")
  dir.create(newIJ)

  file.copy(IJ.path.main, newIJ, recursive=TRUE)
}

identification_list <- paste(folderlist, "/Parallelised_analysis_script_FM_27July2018_Identification.R", sep="")

# UNIX

 ######################################################################
# REFORMAT VIDEOS IF NOT CXD OR AVI
# this has only been tested under Linux
for (foldername in folderlist){

  # if videos are TIFF stacks they need to be converted to avi before bemovi can analyse them
  # the TIFF stacks should be saved in "1_raw_tiff"
  # note: works only for 25 fps
  if(video.format == "tiff"){
    system("mkdir 1_raw")
    system(paste("~/bin/ImageJ/jre/bin/java -Xmx",memory.alloc,"m -jar ~/bin/ImageJ/ij.jar -batch tif_to_avi.ijm",sep=""))
  }
  
  # if videos are .mov videos they need to be converted to avi before bemovi can analyse them
  # the mov videos should be saved in "1_raw_mov"
  if(video.format == "mov"){
    system("mkdir 1_raw")
    # convert all files in the directory
    for (i in 1:length(list.files("1_raw_mov"))){
      # in older distros use ffmpeg with the same syntax
      system(paste("avconv -i 1_raw_mov/",list.files("1_raw_mov")[i]," -f avi -vcodec mjpeg -t ",total_frames/fps," 1_raw/",gsub(".mov", '', list.files("1_raw_mov")[i], ignore.case = T),".avi",sep=""))
    }
  }
}
######################################################################
# TESTING

# check file format and naming
#check_video_file_names(to.data,raw.video.folder,video.description.folder,video.description.file)

# check whether the thresholds make sense (set "dark backgroud" and "red")
#check_threshold_values(to.data, raw.video.folder, ijmacs.folder, 0, difference.lag, thresholds, IJ.path, memory.alloc)

######################################################################
# VIDEO ANALYSIS
#Define number of parallel processes
processes.identifier <- min(detectCores()-1,machine.ram %/% mem.per.identifier)
processes.linker <- min(detectCores()-1,machine.ram %/% mem.per.linker)
processes.overlay <- min(detectCores()-1,machine.ram %/% mem.per.overlay)
folderlist <- paste(folderlist, "/", sep="")

#Function to rscript identification files
script_ident <- function(filename){
  command <- paste("Rscript ", filename, sep="")
  system(command)
}
# identify particles
cl <- makeCluster(processes.identifier)
registerDoParallel(cl)
mclapply(identification_list, script_ident, mc.cores = processes.identifier)
stopCluster(cl)

# link the particles
for (foldername in folderlist){
  to.data <- foldername
  link_particles(to.data, particle.data.folder, trajectory.data.folder, linkrange = trajectory_link_range, disp = trajectory_displacement, start_vid = 1, memory = machine.ram, memory_per_linkerProcess = mem.per.linker)
}

# Do the rest of the steps
for (foldername in folderlist){
  setwd(foldername)
  to.data <- foldername
  command <- paste("mkdir ", foldername, "ijmacs", sep="")
  system(command)
  IJ.path=paste(foldername, "IJ/ImageJ", sep="")
  # merge info from description file and data
  merge_data(to.data, particle.data.folder, trajectory.data.folder, video.description.folder, video.description.file, merged.data.folder)
  
  # load the merged data
  load(paste0(to.data, merged.data.folder, "Master.RData"))
  
  # filter data: minimum net displacement, their duration, the detection frequency and the median step length
  trajectory.data.filtered <- filter_data(trajectory.data, filter_min_net_disp, filter_min_duration, filter_detection_freq, filter_median_step_length)
  
  # summarize trajectory data to individual-based data
  morph_mvt <- summarize_trajectories(trajectory.data.filtered, calculate.median=F, write = T, to.data, merged.data.folder)
  
  # get sample level info
  summarize_populations(trajectory.data.filtered, morph_mvt, write=T, to.data, merged.data.folder, video.description.folder, video.description.file, total_frames)
  
  # create overlays for validation
  create_overlays(trajectory.data.filtered, to.data, merged.data.folder, raw.video.folder, temp.overlay.folder, 
                  overlay.folder, 2048, 2048, difference.lag, type = "label", predict_spec = F, IJ.path=IJ.path, contrast.enhancement = 1, memory = mem.per.overlay)
}

########################################################################
# some cleaning up
for (foldername in folderlist){
  #system("rm -r 2_particle_data")
  #system("rm -r 3_trajectory_data")
  unlink(paste(foldername, "/4a_temp_overlays", sep=""), recursive=TRUE)
  unlink(paste(foldername, "/ijmacs", sep=""), recursive=TRUE)
  unlink(paste(foldername, "/IJ", sep=""), recursive=TRUE)
  ########################################################################
}
