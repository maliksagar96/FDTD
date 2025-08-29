#include <publish.h>
#include <iostream>

Publisher::Publisher() : finished(false) {}

Publisher::~Publisher() {
  stop();
}

void Publisher::start() {
  writer_thread = std::thread(&Publisher::writer_thread_func, this);
}

void Publisher::stop() {
  {
    std::lock_guard<std::mutex> lock(queue_mutex);
    finished = true;
  }
  cv.notify_one();
  if (writer_thread.joinable()) {
    writer_thread.join();
  }
}

void Publisher::push(FieldFrame frame) {
  {
    std::lock_guard<std::mutex> lock(queue_mutex);
    write_queue.push(std::move(frame));
  }
  cv.notify_one();
}

void Publisher::writer_thread_func() {
  while (!finished || !write_queue.empty()) {
    std::unique_lock<std::mutex> lock(queue_mutex);
    cv.wait(lock, [this] { return !write_queue.empty() || finished; });

    while (!write_queue.empty()) {
      auto frame = write_queue.front();
      write_queue.pop();
      lock.unlock();

      for (auto& [name, data] : frame.fields) {
        std::ofstream field_file("data/" + name + "/" + name + std::to_string(frame.t) + ".txt");

        for (int i = 0; i < frame.Nx; ++i) {
          for (int j = 0; j < frame.Ny; ++j) {
            int idx = i * frame.Ny + j;
            field_file << data[idx] << " ";
          }
          field_file << "\n";
        }

        field_file.close();
      }

      lock.lock();
    }
  }
}
