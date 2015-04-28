CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

#the following variables need to be set before invoking the
#dealii autopilot
SET(TARGET "main")
SET(TARGET_SRC 
  burger_settings.cc burger.cc 
  ale_settings.cc ale.cc 
  main.cc
  )

DEAL_II_INVOKE_AUTOPILOT()
