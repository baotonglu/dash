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

#define PROFILE 1

struct System
{
    static void profile(const std::string& name,std::function<void()> body) {
        std::string filename = name.find(".data") == std::string::npos ? (name + ".data") : name;

        // Launch profiler
        pid_t pid;
#ifdef PROFILE
        std::stringstream s;
        s << getpid();
#endif
        int ppid = getpid();
        pid = fork();
        if (pid == 0) {
            // perf to generate the record file
#ifdef PROFILE
            auto fd=open("/dev/null",O_RDWR);
            dup2(fd,1);
            dup2(fd,2);
            exit(execl("/usr/bin/perf","perf","record","-o",filename.c_str(),"-p",s.str().c_str(),nullptr));
#else
            // perf the cache misses of the file
            char buf[200];
            //sprintf(buf, "perf stat -e cache-misses,cache-references,L1-dcache-load-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,r412e -p %d > %s 2>&1",ppid,filename.c_str());
            sprintf(buf, "perf stat -p %d > %s 2>&1",ppid,filename.c_str());
            execl("/bin/sh", "sh", "-c", buf, NULL);
#endif
        }
#ifndef PROFILE
        setpgid(pid, 0);
#endif
        sleep(3);
        // Run body
        body();
        // Kill profiler  
#ifdef PROFILE
        kill(pid, SIGINT);
#else
        kill(-pid, SIGINT);
#endif
        sleep(1);
        //waitpid(pid,nullptr,0);
    }

    static void profile(std::function<void()> body) {
        profile("perf.data",body);
    }
};
