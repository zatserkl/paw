// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "queue.h"
#include "elements.h"
#include "sidelist.h"

//
// class SideList
//

Side* SideList::operator [](int n) const
{
  // Treats list as a closed chain

  if (count == 0) return NULL;   // to prevent cycling
  while (n < 0)      n += count;
  while (n >= count) n -= count;

  Node<Side> *curr = head;
  for (int i=0; i<n; i++) {
    curr=curr->Next;
    if (curr == NULL) curr=head;
  }
  return curr->Item;
}

SideList& SideList::operator =(const SideList& sides)
{
  if (this != &sides) {
    Node<Side> *curr;
    // delete existing list
    while (head != NULL) {
      curr = head;
      head = curr->Next;
      // I have to check this Item in the other nodes in the list!
      delete curr;
      count--;
    }
    // add in list copies of 'queue' items
    for (curr=sides.head; curr; curr=curr->Next) {
      // create a copy of Item
      //
      // I need in typeid stament to create a real copy of Item!
      //
      Side* copy = new Side(*curr->Item);
      Add(copy);
    }
  }
  return *this;
}

void SideList::Reverse()
{
  // In contrast to Queue<T>::Invert head is don't changed,
  // only inverts direction of go round with changing of side ends.

  Node<Side> *tail=head;
  while (tail->Next) tail=tail->Next;
  // connect chain
  tail->Next = head;
  head->Prev = tail;

  Node<Side> *curr=head;
  do {
    // invert the side's ends
    curr->Item->Invert();

    Node<Side> *next = curr->Next;
    curr->Next = curr->Prev;
    curr->Prev = next;
    curr = next;
  } while (curr != head);
  // break chain
  head->Prev->Next = NULL;
  head->Prev = NULL;
}

int SideList::SetHead(Node<Side> *NewHead)
{
  // connect tail with head
  Node<Side> *tail;
  for(tail=head; tail->Next; tail=tail->Next);
  tail->Next = head;
  head->Prev = tail;

  head = NewHead;
  // break chain
  head->Prev->Next = NULL;
  head->Prev = NULL;
  return 0;
}

int SideList::SetHead(Side *side)
{
  Node<Side> *NewHead;
  for (NewHead=head; NewHead; NewHead=NewHead->Next)
    if (NewHead->Item == side) goto FOUND;
  cerr << "ERROR SideList::SetHead(Side *side): "
       << "side did not found!" << endl;
  return 1;

  FOUND:
  // connect tail with head
  Node<Side> *tail;
  for(tail=head; tail->Next; tail=tail->Next);
  tail->Next = head;
  head->Prev = tail;

  head = NewHead;
  // break chain
  head->Prev->Next = NULL;
  head->Prev = NULL;
  return 0;
}

int SideList::SetHead(int n)
{
  if (n == 0) return 0;
  if (n == count) return 0; // in the case of removing of tail side

  // connect tail with head
  Node<Side> *tail;
  for(tail=head; tail->Next; tail=tail->Next);
  tail->Next = head;
  head->Prev = tail;

  Node<Side> *NewHead=head;
  for (int i=0; i<n; i++) NewHead=NewHead->Next;
  head = NewHead;
  // break chain
  head->Prev->Next = NULL;
  head->Prev = NULL;
  return 0;
}

int SideList::SetHeadAfter(Side *side)
{
  Node<Side> *SideNode;
  for (SideNode=head; SideNode; SideNode=SideNode->Next)
    if (SideNode->Item == side) goto FOUND;
  cerr << "ERROR SideList::SetHeadAfter(Side *side): "
       << "side did not found!" << endl;
  return 1;

  FOUND:
  // find tail of list
  Node<Side> *tail;
  for (tail=head; tail->Next; tail=tail->Next);
  if (SideNode != tail)
  {
    // connect tail with head
    tail->Next = head;
    head->Prev = tail;

    head = SideNode->Next;
    // break chain
    head->Prev->Next = NULL;
    head->Prev = NULL;
  }
  return 0;
}

int SideList::SetHeadAfter(int n)
{
  return SetHeadAfter(GetItem(n));
}

int SideList::SetHeadBefore(Side *side)
{
  Node<Side> *SideNode;
  for (SideNode=head; SideNode; SideNode=SideNode->Next)
    if (SideNode->Item == side) goto FOUND;
  cerr << "ERROR SideList::SetHeadBefore(Side *side): "
       << "side did not found!" << endl;
  return 1;

  FOUND:
  // connect tail with head
  Node<Side> *tail;
  for(tail=head; tail->Next; tail=tail->Next);
  tail->Next = head;
  head->Prev = tail;

  head = SideNode->Prev;
  // break chain
  head->Prev->Next = NULL;
  head->Prev = NULL;
  return 0;
}

int SideList::SetHeadBefore(int n)
{
  return SetHeadBefore(GetItem(n));
}

void SideList::Move(double dx, double dy, double dtheta, int inv)
{
  for (Node<Side> *curr = head; curr; curr=curr->Next) {
    curr->Item->Move(dx,dy,dtheta,inv);
  }
}

void SideList::Move(const Point& dr)
{
  for (Node<Side> *curr = head; curr; curr=curr->Next) {
    curr->Item->Move(dr);
  }
}

void SideList::Optimize()
{
  // Perform the finish optimizing of the side list

  if (count <= 1) return;

  Node<Side> *tail=head;
  while (tail->Next) tail=tail->Next;

  // Eliminate the done sequences

  Node<Side> *curr;
  for (curr=head; curr; curr=curr->Next) {
    if (curr->Item->status != todo)
    {
      Node<Side> *done1,*done2,*donecurr;
      done1 = curr;
      donecurr = done1;
      while (donecurr != NULL && donecurr->Item->status != todo) {
        done2 = donecurr;
        donecurr = donecurr->Next;
      }
      if (done2 == tail) {
        // mark tail tile(s) as an 'extra'
        for (curr=done1; curr != done2->Next; curr=curr->Next)
        {
          curr->Item->status = extra;
        }
        break; // we have reached the tail, now curr==NULL
      }
      else if (done2 != done1)
      {
        // we have two or more done sides, mark they 'extra'
        for (curr=done1; curr != done2->Next; curr=curr->Next)
        {
          curr->Item->status = extra;
        }
        // insert a 'lift' side
        Side *liftside = new Side(done1->Item->p1,done2->Item->p2, lift);
        InsertBefore(done1->Item, liftside);
      }
    }
  }
}

void SideList::Show() const
{
  for (Node<Side> *curr=head; curr; curr=curr->Next) {
    curr->Item->Show();
    cout << endl;
  }
}

void SideList::Rotate(const Point& O, double theta) {
  for (Node<Side> *curr=head; curr; curr=curr->Next) {
    curr->Item->Rotate(O,theta);
  }
}

void SideList::TurnBack(double A, double B) {
  for (Node<Side> *curr=head; curr; curr=curr->Next) {
    curr->Item->TurnBack(A,B);
  }
}
