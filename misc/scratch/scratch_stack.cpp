  // stack size
  //
  // attempt to avoid seg fault on lsu::wru3 startup on Cygwin
  // http://stackoverflow.com/questions/2275550/change-stack-size-for-a-c-application-in-linux-during-compilation-with-gnu-com
  {
    const rlim_t kStackSize = 64 * 1024 * 1024;   // min stack size = 16 MB
    struct rlimit rl;
    int result;
    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
      {
        std::cout << "initial stack size (Mb) " << rl.rlim_cur/double(1024*1024) << std::endl;
        if (rl.rlim_cur < kStackSize)
          {
            rl.rlim_cur = kStackSize;
            std::cout << "resetting stack size (Mb) " << rl.rlim_cur/double(1024*1024) << std::endl;
            result = setrlimit(RLIMIT_STACK, &rl);
            if (result != 0)
              {
                fprintf(stderr, "setrlimit returned result = %d\n", result);
              }
          }
      }
    result = getrlimit(RLIMIT_STACK, &rl);
    if (result == 0)
      {
        std::cout << "current stack size (Mb) " << rl.rlim_cur/double(1024*1024) << std::endl;
      }  
  }
