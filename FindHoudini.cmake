# Find Houdini
if (MSVC)
    set(Houdini_PATH "C:/Program Files/Side Effects Software/Houdini 20.0.697")
elseif (APPLE)
    set(Houdini_PATH "/Applications/Houdini/Houdini20.0.697/Frameworks/Houdini.framework/Versions/20.0/Resources")
endif ()
set(Houdini_DIR ${Houdini_PATH}/toolkit/cmake)
find_package(Houdini REQUIRED)
