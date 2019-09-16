#pragma once

#include <sstream>
#include <iostream>
#include <functional>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

struct System
{
    static void profile(const std::string& name,std::function<void()> body) {
        std::string filename = name.find(".data") == std::string::npos ? (name + ".data") : name;

        // Launch profiler
        pid_t pid;
        //std::stringstream s;
        //s << getpid();
        int ppid = getpid();
        pid = fork();
        if (pid == 0) {
            /*
            // perf to generate the record file
            auto fd=open("/dev/null",O_RDWR);
            dup2(fd,1);
            dup2(fd,2);
            exit(execl("/usr/bin/perf","perf","record","-o",filename.c_str(),"-p",s.str().c_str(),nullptr));
            */
            // perf the cache misses of the file
            char buf[200];
            sprintf(buf, "perf stat -e cache-misses,cache-references,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,r412e -p %d > %s 2>&1",ppid,filename.c_str());
            execl("/bin/sh", "sh", "-c", buf, NULL);
        }
        setpgid(pid, 0);
        sleep(1);
        // Run body
        body();

        // Kill profiler  
        kill(-pid,SIGINT);
        sleep(1);
        //waitpid(pid,nullptr,0);
    }

    static void profile(std::function<void()> body) {
        profile("perf.data",body);
    }
};