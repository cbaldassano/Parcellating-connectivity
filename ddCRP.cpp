#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <iostream>

using namespace std;

void ThreadFunction(int i)
{
  boost::this_thread::sleep(boost::posix_time::milliseconds(i*1000));
  cout << i << endl;
  if (i < 10) {
    boost::thread_group Tg;
    Tg.create_thread(boost::bind(ThreadFunction,i+10));
    Tg.join_all();
  }
}

int main()
{
  boost::thread_group Tg;

  for (int i = 0; i < 5; ++i)
    Tg.create_thread(boost::bind(ThreadFunction, i));

  Tg.join_all();
  cout << "Done!" << endl;
  return 0;
}
