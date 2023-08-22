  include(FetchContent)
  FetchContent_Declare(
      libPerm
      GIT_REPOSITORY      https://github.com/Krzmbrzl/libPerm.git
      GIT_TAG             v1.5.0
  )
  FetchContent_MakeAvailable(libPerm)

