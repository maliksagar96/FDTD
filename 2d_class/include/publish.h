#pragma once
#include <map>
#include <string>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <thread> 

// A single snapshot of fields at a given timestep
struct FieldFrame {
  int t;
  int Nx, Ny;
  std::map<std::string, std::vector<double>> fields;
};

class Publisher {
public:
  Publisher();
  ~Publisher();

  void start();
  void stop();
  void push(FieldFrame frame);

private:
  void writer_thread_func();
  std::queue<FieldFrame> write_queue;
  std::mutex queue_mutex;
  std::condition_variable cv;
  bool finished;
  std::thread writer_thread;
};
