######################################################################
# R script for analysing video files with BEMOVI (www.bemovi.info)
#
# Felix Moerman
#
# July 2018
######################################################################
rm(list=ls())

# load package
#library(devtools)
#install_github("efronhofer/bemovi", ref="experimental")
library(bemovi)

# UNIX
# set paths to ImageJ and particle linker standalone
IJ.path.main <- "IJ/ImageJ"
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
mem.per.identifier <- c(30000)

to.data <- paste(getwd(),"/",sep="")
IJ.path.main <- paste(to.data, "IJ/ImageJ", sep="")
# identify particles
locate_and_measure_particles(to.data, raw.video.folder, particle.data.folder, difference.lag, 
         thresholds, min_size = particle_min_size, max_size = particle_max_size, 
         IJ.path=IJ.path.main, memory=mem.per.identifier)
