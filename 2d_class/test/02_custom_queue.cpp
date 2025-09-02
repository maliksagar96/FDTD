#include <iostream>


using namespace std;

class Node{
  public:
    int val;
    Node *next;

    Node():val(0), next(nullptr) {}
    Node(int x):val(x), next(nullptr) {}
    Node(int x, Node *next):val(x), next(next){}
};

class queue {
  
  Node *front, *rear;
  int size;

  public:

    void push() {
      if(!front) 
      if(!rear) 
    }

    int pop() {
      return -1;
    }

    int max() {
      return -1;
    }

    int min() {
      return -1;
    }

    bool isEmpty() {
      return size == 0;
    }

    int getSize() {
      return size;
    }
    
};




int main() {


  return 0;
}
