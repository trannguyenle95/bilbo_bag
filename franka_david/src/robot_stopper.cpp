#include <ros/ros.h>
#include <std_msgs/Bool.h>
#include <fstream>
#include <stdio.h>
#include <typeinfo>

#include <control_remote.h>

#include <chrono>


//NOTE: file-based UDP stopper

//from: https://www.geeksforgeeks.org/udp-server-client-implementation-c/
 
// Server side implementation of UDP client-server model
#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <chrono>
 
#define PORT     8080
#define MAXLINE 1 //flag just one char

// Driver code
int main() {

nice(-20);

//file-based stopper
bool franka3 = true; //set which robot the script is running on
std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
std::string outfile = HOME + "/catkin_ws/src/Data/fr3_error.txt";
//std:string topic = "/franka3_control_node/franka3_state"; //franka2 subscribes to error franka3_state
if (franka3){
    outfile = HOME + "/catkin_ws/src/Data/fr2_error.txt";
    //topic = "/franka2_control_node/franka2_state"; //franka3 subscribes to error franka2_state
}
remove(outfile.c_str());


int sockfd;
char buffer[MAXLINE];
struct sockaddr_in servaddr, cliaddr;
 
// Creating socket file descriptor
if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) {
perror("socket creation failed");
exit(EXIT_FAILURE);
}
 
memset(&servaddr, 0, sizeof(servaddr));
memset(&cliaddr, 0, sizeof(cliaddr));
 
// Filling server information
servaddr.sin_family    = AF_INET; // IPv4
//servaddr.sin_addr.s_addr = INADDR_ANY;
servaddr.sin_port = htons(PORT);
servaddr.sin_addr.s_addr = inet_addr("130.233.123.182");
 
// Bind the socket with the server address
if ( bind(sockfd, (const struct sockaddr *)&servaddr,  
sizeof(servaddr)) < 0 )
{
perror("bind failed");
exit(EXIT_FAILURE);
}
 
socklen_t len;
int n;
 
len = sizeof(cliaddr);  //len is value/result
 
// n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
// MSG_WAITALL, ( struct sockaddr *) &cliaddr,
// &len);



//from: https://stackoverflow.com/questions/2876024/linux-is-there-a-read-or-recv-from-socket-with-timeout
struct timeval tv;
tv.tv_sec = 0; //timeout_in_seconds
tv.tv_usec = 2000; //20 000 microseconds = 20 ms, 10ms = stops at 14-16 ms, 2ms = stops after 4 ms - without stop button // 5ms = stops after 8ms - without stop button // 8 ms = 12ms

n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
0, ( struct sockaddr *) &cliaddr,
&len);
std::cout << "message received" << n << std::endl;
std::cout << "buffer[0] " << buffer[0] << std::endl;
bool test = (buffer[0] == '0');
std::cout << "buffer[0] == '0'" << test  << std::endl;
setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, (const char*)&tv, sizeof tv);

bool reset = 1;

int i = 0;

while(1){
  n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
  0, ( struct sockaddr *) &cliaddr,
  &len);
  //std::cout << "buffer[0] " << buffer[0] << std::endl;

  //string flag = buffer[0];

  if(buffer[0]=='1'){
    std::cout << "reset received" << n << std::endl;
    remove(outfile.c_str());
    reset = 1;
  }

  if((n <= 0) && (reset == 1)){
    std::ofstream output(outfile); //ADDED for file-based stopper
    std::cout << "stopped" << n << std::endl;
    reset = 0;
  }

  //buffer[n] = '\0';
  //std::cout << "message received " << buffer << std::endl;
  i++;
  if((i%5000)==0){
    std::cout << "i " << i << std::endl;
  }
}





while(n>0){
  //std::cout << std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
  n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
  0, ( struct sockaddr *) &cliaddr,
  &len);
  buffer[n] = '\0';
  std::cout << "message received " << buffer << std::endl;
}

//buffer[n] = '\0';
//printf("Client : %s\n", buffer);
//std::cout << "timeout" << n << " time " << std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
std::ofstream output(outfile); //ADDED for file-based stopper

 
return 0;
}



// int main(int argc, char **argv)
// {
//   bool franka3 = true; //set which robot the script is running on

//   ros::init(argc, argv, "listener");
//   //stop node when error message is received
//   std_msgs::Bool msg;
//   std:string topic = "/franka3_control_node/franka3_state"; //franka2 subscribes to error franka3_state
//   if (franka3){
//     topic = "/franka2_control_node/franka2_state"; //franka3 subscribes to error franka2_state
//   }

//   msg  = *(ros::topic::waitForMessage<std_msgs::Bool>(topic));

//   return 0;
// }


//NOTE: file-based stopper
// class Stopper
// {
//   public:
//   //uintptr_t _address;
//   std::string _file;
//   void callback(const std_msgs::Bool::ConstPtr& msg)
//   {
//       std::cout << "callback for stopping" << std::endl;
//       std::cout << std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
//       std::ofstream output(_file);
//   }

// };

// int main(int argc, char **argv)
// {
//   bool franka3 = true; //set which robot the script is running on

//   std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
//   std::string outfile = HOME + "/catkin_ws/src/Data/fr3_error.txt";
//   std:string topic = "/franka3_control_node/franka3_state"; //franka2 subscribes to error franka3_state
//   if (franka3){
//     outfile = HOME + "/catkin_ws/src/Data/fr2_error.txt";
//     topic = "/franka2_control_node/franka2_state"; //franka3 subscribes to error franka2_state
//   }

//   remove(outfile.c_str());
//   ros::init(argc, argv, "listener");

//   Stopper stopper;
//   ros::NodeHandle n("~");
//   stopper._file = outfile;

//   //ros::Subscriber sub = n.subscribe(topic, 1000, &Stopper::callback, &stopper);
//   ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1, &Stopper::callback, &stopper, ros::TransportHints().maxDatagramSize(1).udp());
//   ros::spin();
  
//   return 0;
// }


//NOTE: kill via pid
// #include <sys/types.h>
// #include <signal.h>

// class Stopper
// {
//   public:
//   //uintptr_t _address;
//   int _address;
//   void callback(const std_msgs::Bool::ConstPtr& msg)
//   {
//     kill(_address, SIGTERM);
//   }

// };



// int main(int argc, char **argv)
// {
//   // std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
//   // std::string outfile = HOME + "/catkin_ws/src/Data/fr2_error.txt";
//   // remove(outfile.c_str());
//   ros::init(argc, argv, "listener");

//   Stopper stopper;
//   int address;
//   ros::NodeHandle n("~");
//   n.getParam("address", address);
//   stopper._address = address;
//   std::cout << "received address: " << address << std::endl;
//   ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1000, &Stopper::callback, &stopper);
//   ros::spin();

//   return 0;
// }

//STOP from same PC, check timing
// class Stopper
// {
//   public:
//   //uintptr_t _address;
//   std::string _address;
//   void callback(const std_msgs::Bool::ConstPtr& msg)
//   {
//       //exit(0);
//       std::cout << "callback for stopping" << std::endl;
      
//       if(msg->data){
//           // //Using this: https://stackoverflow.com/questions/47591008/how-do-i-convert-a-string-to-an-address
//           // // std::stringstream ss;
//           // // ss << _address;
//           // // std::uintptr_t temp;
//           // // ss >> std::hex >> temp;
//           // // //std::cout << "address: " << ss << std::endl; 
//           // // std::cout << "here" << std::endl;
//           // // *(reinterpret_cast<bool*>(temp)) = false;

//           // StopFlag flag_object; // = new StopFlag();
//           // // prompt the id of this->_flag_object->_flag;
//           // std::cout << "flag value: " << flag_object._flag << std::endl;
//           // std::cout << "flag address: " << &(flag_object._flag) << std::endl;
//           // //std::string ss = &(flag_object->_flag);


//           // bool qs = false;
//           // std::stringstream ss;
//           // ss << &(flag_object._flag);
//           // //ss << &qs;
//           // //std::cout << "qs type: " << typeid(&qs).name() << std::endl;
//           // std::cout << "ss: " << ss.str() << std::endl;
//           // std::stringstream ss2;
//           // ss2 << _address;
//           // //std::cout << "ss2 type: " << typeid(_address).name() << std::endl;
//           // std::cout << "ss2: " << ss2.str() << std::endl;

//           // std::uintptr_t temp;
//           // ss2 >> std::hex >> temp;
//           // bool newValue = true;
//           // std::cout << "here" << std::endl;
//           // //std::cout << "current value:" << *(reinterpret_cast<bool*>(temp)) << std::endl;
//           // //*(reinterpret_cast<bool*>(temp)) = newValue;
//           // StopFlag* pp = reinterpret_cast<StopFlag*>(temp);
//           // pp->_flag = newValue;
          
//           // //int* pp = reinterpret_cast<int*>(temp);

//           std::cout << "stop flag received" << std::endl;
//           std::cout << std::chrono::duration_cast< std::chrono::milliseconds >(std::chrono::system_clock::now().time_since_epoch()).count() << std::endl;
          

//       }
//   }

// };



// int main(int argc, char **argv)
// {
//   // std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
//   // std::string outfile = HOME + "/catkin_ws/src/Data/fr2_error.txt";
//   // remove(outfile.c_str());
//   ros::init(argc, argv, "listener");
//   //ros::Rate rate(1000);

//   Stopper stopper;
//   std::string address = "default";
//   ros::NodeHandle n("~");
//   n.getParam("address", address);
//   stopper._address = address;
//   std::cout << "received address: " << address << std::endl;
//   //ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1, &Stopper::callback, &stopper);
//   //ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1, &Stopper::callback, &stopper, ros::TransportHints().reliable().tcpNoDelay());
//   //ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1, &Stopper::callback, &stopper, ros::TransportHints().tcpNoDelay());
//   ros::Subscriber sub = n.subscribe("/franka2_control_node/franka2_state", 1, &Stopper::callback, &stopper, ros::TransportHints().maxDatagramSize(1).udp());
//   ros::spin();


//   //stop node when error message is received
//   // std_msgs::Bool msg;
//   // msg  = *(ros::topic::waitForMessage<std_msgs::Bool>("/franka2_control_node/franka2_state"));

//   return 0;
// }