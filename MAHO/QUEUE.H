// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#include <stdlib.h>
#include <iostream.h>

template <class T>
struct Node {   // the list item
  T* Item;      // can be T or any class derived from T
  Node *Next;   // point to next Node object
  Node *Prev;   // point to prev Node object

  Node(T* NewItem, Node<T>* prev, Node<T>* next)
                   {Item=NewItem; Prev=prev; Next=next;}
  Node(T* NewItem) {Item=NewItem; Prev=NULL; Next=NULL;}
  Node()           {Item=NULL;    Prev=NULL; Next=NULL;}
};

template <class T>
class Queue {

// Only one instance of item can be stored

protected:
  int count;
  Node<T>* head;
  Node<T>* GetNext(Node<T>* node) const;
public:
  Queue();
  //Queue(const Queue<T>& queue);
  ~Queue();
  virtual T* operator [](int n) const;
  Queue<T>&  operator =(const Queue<T>& queue);
  virtual Queue<T>& Add(T* NewItem);
  virtual Queue<T>& operator << (T* NewItem); // without 'virtual' linker error???
  int Count() const {return count;}
  int Index(T* Item) const;
  int Exist(T* Item) const;
  T* GetItem(int n) const;
  T* Head() const {return head->Item;}
  T* Tail() const;
  int InsertAfter(T* TheItem,  T* NewItem);
  int InsertAfter(int n,       T* NewItem);
  int InsertBefore(T* TheItem, T* NewItem);
  int InsertBefore(int n,      T* NewItem);
  T* Remove(T* TheItem);
  T* Remove(int n);
  int Kill(Node<T> *TheNode);
  int Kill(T* TheItem);
  int Kill(int n);
  int ReplaceItem(T *Item, T *ByItem);
  int ReplaceItem(int n, T *ByItem);
  int Swap(T* Item1, T* Item2);
  int Swap(int n1, int n2);
  void Invert();
  //void Applay(void (T::*pmf)());
  void Applay(void (*pf)(T& member));
};

template <class T> Queue<T>::Queue()
{
  count = 0;
  head = NULL;
}

/*
template <class T> Queue<T>::Queue(const Queue<T>& queue)
{
  for(Node<T> *curr=queue.head; curr; curr=curr->Next) Add(curr->Item);
}
*/

template <class T> Queue<T>::~Queue()  // destructor
{
  while (head != NULL) {
    Node<T>* curr = head;
    head = curr->Next;
    delete curr->Item;
    delete curr;
    count--;
  }
}

//
// very dangerous operator without typeid statment
//
template <class T>
Queue<T>& Queue<T>::operator =(const Queue<T>& queue)
{
  const char routine[] = "Queue<T>::operator =(const Queue<T>& queue)";

  cout << routine << endl;

  if (this != &queue) {
    Node<T> *curr;
    // delete existing list
    while (head != NULL) {
      curr = head;
      head = curr->Next;
      delete curr;
      count--;
    }
    // add in list copies of 'queue' items
    for (curr=queue.head; curr; curr=curr->Next) {
      // create a copy of Item
      //
      // I need in typeid statment to create a real copy of Item!
      //

      T* copy = new T(*curr->Item);
      Add(copy);
      //Add(curr->Item);
    }
  }
  return *this;
}

template <class T>
T* Queue<T>::operator[](int n) const
{
  const char routine[] = "Queue<T>::operator[](int n)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!";
    exit(EXIT_FAILURE);
  }
  Node<T>* curr = head;
  for (int i=0; i<n; i++) curr=curr->Next;
  return curr->Item;
}

template <class T>
T* Queue<T>::Tail() const
{
  Node<T>* tail=head;
  while (tail->Next) tail=tail->Next;
  return tail->Item;
}

template <class T>
T* Queue<T>::GetItem(int n) const
{
  const char routine[] = "Queue<T>::GetItem(int n)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!";
    exit(EXIT_FAILURE);
  }
  Node<T>* curr = head;
  for (int i=0; i<n; i++) curr=curr->Next;
  return curr->Item;
}

template <class T>
int Queue<T>::Index(T* Item) const
{
  const char routine[] = "Queue<T>::Index(T* Item)";

  int n=0;
  for (Node<T>* curr=head; curr; curr=curr->Next) {
    if (curr->Item == Item) goto FOUND;
    n++;
  }
  cerr << routine << ": Item did not found!"
       << endl;
  return -1;

  FOUND:
  return n;
}

template <class T>
int Queue<T>::Exist(T* Item) const
{
  for (Node<T>* curr=head; curr; curr=curr->Next) {
    if (curr->Item == Item) return 1;
  }
  return 0;
}

template <class T>
T* Queue<T>::Remove(T* TheItem)
{
  // return TheItem in success or NULL then TheItem did not found

  //cout << routine << endl;

  const char routine[] = "Queue<T>::Remove(T* TheItem)";

  Node<T>* curr=head;
  while (curr) {
    if (curr->Item == TheItem) goto FOUND;
    curr=curr->Next;
  }
  cerr << routine << ": TheItem did not found!"
       << endl;
  return NULL;

  FOUND:
  if (curr != head) {
    curr->Prev->Next = curr->Next;
    // if curr isn't a tail
    if (curr->Next != NULL) curr->Next->Prev = curr->Prev;
  }
  else {
    head = head->Next;
    if (head) head->Prev = NULL;
  }

  delete curr;
  count--;
  return TheItem;
}

template <class T>
T* Queue<T>::Remove(int n)
{
  // return pointer to item in success or NULL then n is out of range

  //cout << "Queue<T>::Remove(int n)" << endl;

  const char routine[] = "Queue<T>::Remove(int n)";

  if (n >= count || n < 0) {
    cerr << routine << ": n is out of range!"
         << endl;
    return NULL;
  }

  Node<T>* curr=head;
  for (int i=0; i<n; i++) curr=curr->Next;
  T* TheItem = curr->Item;

  if (curr != head) {
    curr->Prev->Next = curr->Next;
    // if curr isn't a tail
    if (curr->Next != NULL) curr->Next->Prev = curr->Prev;
  }
  else {
    head = head->Next;
    if (head) head->Prev = NULL;
  }

  delete curr;
  count--;
  return TheItem;
}

template <class T>
int Queue<T>::Kill(Node<T> *TheNode)
{
  // return 1 in success or 0 then TheNode did not found or queue is empty

  const char routine[] = "Queue<T>::Kill(Node<T> *TheNode)";

  Node<T>* curr=head;
  while (curr) {
    if (curr == TheNode) goto FOUND;
    curr=curr->Next;
  }
  cerr << routine << ": TheNode did not found!"
       << endl;
  return 0;

  FOUND:;
  if (curr != head) {
    curr->Prev->Next = curr->Next;
    // if curr isn't a tail
    if (curr->Next != NULL) curr->Next->Prev = curr->Prev;
  }
  else {
    head = head->Next;
    if (head) head->Prev = NULL;
  }

  // Kill item
  delete TheNode->Item;
  delete TheNode;
  count--;
  return 1;
}

template <class T>
int Queue<T>::Kill(T* TheItem)
{
  // return 1 in success or 0 then TheItem did not found or queue is empty

  const char routine[] = "Queue<T>::Kill(T* TheItem)";

  Node<T>* curr=head;
  while (curr) {
    if (curr->Item == TheItem) goto FOUND;
    curr=curr->Next;
  }
  cerr << routine << ": TheItem did not found!"
       << endl;
  return 0;

  FOUND:;
  if (curr != head) {
    curr->Prev->Next = curr->Next;
    // if curr isn't a tail
    if (curr->Next != NULL) curr->Next->Prev = curr->Prev;
  }
  else {
    head = head->Next;
    if (head) head->Prev = NULL;
  }

  delete curr->Item;
  delete curr;
  count--;
  return 1;
}

template <class T>
int Queue<T>::Kill(int n)
{
  // return 1 in success or 0 then n is out of range or queue is empty

  const char routine[] = "Queue<T>::Kill(int n)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!";
    exit(EXIT_FAILURE);
  }

  Node<T>* curr=head;
  for (int i=0; i<n; i++) curr=curr->Next;

  if (curr != head) {
    curr->Prev->Next = curr->Next;
    // if curr isn't a tail
    if (curr->Next != NULL) curr->Next->Prev = curr->Prev;
  }
  else {
    head = head->Next;
    if (head) head->Prev = NULL;
  }

  delete curr->Item;
  delete curr;
  count--;
  return 0;
}

template <class T>
Queue<T>& Queue<T>::Add(T* NewItem)
{
  const char routine[] = "Queue<T>::Add(T* NewItem)";

  if (head == NULL) head = new Node<T>(NewItem,NULL,NULL);
  else {
    // add at the end
    Node<T>* tail = head;
    while (tail->Next != NULL) {
      if (tail->Item == NewItem) {
        cerr << routine << ": NewItem already stored" << endl;
        exit(1);
      }
      tail=tail->Next;
    }
    tail->Next = new Node<T>(NewItem,tail,NULL);
  }
  count++;
  return *this;
}

template <class T>
Queue<T>& Queue<T>:: operator << (T* NewItem)
{
  const char routine[] = "Queue<T>:: operator << (T* NewItem)";

  // Add(NewItem);
  if (head == NULL) head = new Node<T>(NewItem,NULL,NULL);
  else {
    Node<T>* tail = head;
    while (tail->Next) {
      if (tail->Item == NewItem) {
        cerr << routine << " : NewItem already stored" << endl;
        exit(1);
      }
      tail=tail->Next;
    }
    tail->Next = new Node<T>(NewItem,tail,NULL);
  }
  count++;
  return *this;
}

template <class T>
int Queue<T>::InsertAfter(T* TheItem, T* NewItem)
{
  // return Index of inserted item
  // or -1 then TheItem did not found

  const char routine[] = "Queue<T>::InsertAfter(T* TheItem, T* NewItem)";

  int i=0;
  Node<T>* curr=head;
  while (curr) {
    if (curr->Item == TheItem) goto FOUND;
    curr=curr->Next;
    i++;
  }
  cerr << routine << ": TheItem did not found!"
       << endl;
  return -1;

  FOUND:;
  if (Exist(NewItem)) {
    cerr << routine << ": NewItem already stored" << endl;
    exit(1);
  }
  Node<T>* N = new Node<T>(NewItem, curr, curr->Next);
  curr->Next = N;
  if (N->Next) N->Next->Prev = N;
  i++;
  count++;
  return i;
}

template <class T>
int Queue<T>::InsertAfter(int n, T* NewItem)
{
  // return Index of inserted item or -1 then n is out of range
  // Note: n<0 means insert before head

  const char routine[] = "Queue<T>::InsertAfter(int n, T* NewItem)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!" << endl;
    return -1;
  }

  if (Exist(NewItem)) {
    cerr << routine << ": NewItem already stored" << endl;
    exit(1);
  }

  int i;
  Node<T>* curr=head;
  for (i=0; i<n; i++) curr=curr->Next;

  Node<T>* N = new Node<T>(NewItem, curr, curr->Next);
  curr->Next = N;
  if (N->Next) N->Next->Prev = N;
  i++;
  count++;
  return i;
}

template <class T>
int Queue<T>::InsertBefore(T* TheItem, T* NewItem)
{
  // return Index of inserted item
  // or -1 then TheItem did not found

  const char routine[] = "Queue<T>::InsertBefore(T* TheItem, T* NewItem)";

  int i=0;
  Node<T>* curr=head;
  while (curr) {
    if (curr->Item == TheItem) goto FOUND;
    curr=curr->Next;
    i++;
  }
  cerr << routine << ": TheItem did not found!" << endl;
  return -1;

  FOUND:;
  if (Exist(NewItem)) {
    cerr << routine << ": NewItem already stored" << endl;
    exit(1);
  }
  // insert node before curr
  Node<T>* N = new Node<T>(NewItem, curr->Prev, curr);
  if (curr != head) {
    curr->Prev->Next = N;
    // index of inserted item is equal to index of curr;
  }
  else {
     head = N;
     i = 0;
  }
  curr->Prev = N;
  count++;
  return i;
}

template <class T>
int Queue<T>::InsertBefore(int n, T* NewItem)
{
  // return Index of inserted item or -1 then n is out of range
  // Note: n<0 means insert before head

  const char routine[] = "Queue<T>::InsertBefore(int n, T* NewItem)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!" << endl;
    return -1;
  }

  if (Exist(NewItem)) {
    cerr << routine << ": NewItem already stored" << endl;
    exit(1);
  }

  int i;
  Node<T> *curr=head;
  for (i=0; i<n; i++) curr=curr->Next;

  Node<T>* N = new Node<T>(NewItem, curr->Prev, curr);
  if (curr != head) {
    curr->Prev->Next = N;
    // index of inserted item is equal to index of curr;
  }
  else {
     head = N;
     i = 0;
  }
  curr->Prev = N;
  count++;
  return i;
}

template <class T>
int Queue<T>::ReplaceItem(T *Item, T *ByItem)
{
  const char routine[] = "Queue<T>::ReplaceItem(T *Item, T *ByItem)";

  if (Exist(ByItem)) {
    cerr << routine << ": ByItem already stored" << endl;
    exit(1);
  }
  Node<T> *curr;
  for (curr=head; curr; curr=curr->Next) {
    if(curr->Item == Item) {
      curr->Item = ByItem;
      return 1;
    }
  }
  cerr << routine << ": Item did not found!" << endl;
  return 0;
}

template <class T>
int Queue<T>::ReplaceItem(int n, T *ByItem)
{
  const char routine[] = "Queue<T>::ReplaceItem(int n, T *ByItem)";

  if (n < 0 || n >= count) {
    cerr << routine << ": n is out of range!" << endl;
    return 0;
  }

  if (Exist(ByItem)) {
    cerr << routine << ": ByItem already stored" << endl;
    exit(1);
  }

  Node<T> *curr=head;
  for (int i=0; i<n; i++) curr=curr->Next;
  curr->Item = ByItem;
  return 1;
}

template <class T>
int Queue<T>::Swap(T* Item1, T* Item2)
{
  const char routine[] = "Queue<T>::Swap(T* Item1, T* Item2)";

  Node<T> *p1,*p2;
  for (p1=head; p1; p1=p1->Next) if(p1->Item == Item1) goto FOUND1;
  cerr << routine << ": Item1 did not found!" << endl;
  return 0;

  FOUND1:;
  for (p2=head; p2; p2=p2->Next) if(p2->Item == Item2) goto FOUND2;
  cerr << routine << ": Item2 did not found!" << endl;
  return 0;

  FOUND2:;
  p1->Item = Item2;
  p2->Item = Item1;
  return 1;
}

template <class T>
int Queue<T>::Swap(int n1, int n2)
{
  const char routine[] = "Queue<T>::Swap(int n1, int n2)";

  if (n1 < 0 || n1 >= count) {
    cerr << routine << ": n1 is out of range!" << endl;
    return 0;
  }
  if (n2 < 0 || n2 >= count) {
    cerr << routine << ": n2 is out of range!" << endl;
    return 0;
  }

  Node<T> *p1,*p2;
  p1 = p2 = head;
  int i;
  for (i=0; i<n1; i++) p1=p1->Next;
  for (i=0; i<n2; i++) p2=p2->Next;

  T* buf = p1->Item;
  p1->Item = p2->Item;
  p2->Item = buf;
  return 1;
}

template <class T>
void Queue<T>::Invert()
{
  Node<T>* curr=head;
  while (curr) {
    Node<T>* next = curr->Next;
    curr->Next = curr->Prev;
    curr->Prev = next;
    // shift head
    head = curr;
    // pointer to next node now stored in 'next'
    curr = next;
  }
}

/*
template <class T>
void Queue<T>::Applay(void (T::*pmf)())
{
   Node<T>* curr = head;
   while (curr != NULL) {
     (curr->Item->*pmf)();
     curr = curr->Next;
   }
}
*/

template <class T>
void Queue<T>::Applay(void (*pf)(T& member))
{
   Node<T>* curr = head;
   while (curr != NULL) {
     (*pf)(*curr->Item);
     curr = curr->Next;
   }
}

