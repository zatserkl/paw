// Andriy_Zatserklyaniy@fnal.gov

template <class T> class Queue;
struct Point;
struct Side;

class SideList: public Queue<Side> {
//
// Main difference between SideList and Queue<T>:
// for some operations SideList treats list as CLOSED chain
//
protected:
  int SetHead(Node<Side> *NewHead);
public:
  Side* operator [](int n) const;
  SideList& operator =(const SideList& sides);
  void Reverse();
  void Move(double dx, double dy, double dtheta, int inv);
  void Move(const Point& dr);
  void Rotate(const Point& O, double theta);
  void TurnBack(double A, double B);
  int SetHead(Side *side);
  int SetHead(int n);
  int SetHeadAfter(Side *side);
  int SetHeadAfter(int n);
  int SetHeadBefore(Side *side);
  int SetHeadBefore(int n);
  void Optimize();
  void Show() const;
};
