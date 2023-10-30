#include <ros/ros.h>
#include <std_msgs/Bool.h>
#include <fstream>
#include <stdio.h>
#include <typeinfo>
#include <control_remote.h>
#include <chrono>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>


//file-based UDP stopper

//from UDP server based on: https://www.geeksforgeeks.org/udp-server-client-implementation-c/

#define PORT     8080
#define MAXLINE 1 //flag just one char

// Driver code
int main() {

  bool franka3 = true; //set which robot the script is running on

  if(franka3){
    std::cout << "Stopper for Franka3" << std::endl;
  } else {
    std::cout << "Stopper for Franka2" << std::endl;
  }

  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  std::string outfile = HOME + "/catkin_ws/src/Data/fr3_error.txt";
  if (franka3){
      outfile = HOME + "/catkin_ws/src/Data/fr2_error.txt";
  }
  remove(outfile.c_str()); //remove file that flags stopping to allow the robot to run


  int sockfd;
  char buffer[MAXLINE];
  struct sockaddr_in servaddr, cliaddr;
  
  // Creating socket file descriptor
  if ( (sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0 ) {
  perror("socket creation failed");
  exit(EXIT_FAILURE);
  }

  const int toggle = 1;
  setsockopt(sockfd, SOL_SOCKET, SO_REUSEPORT, &toggle, sizeof(toggle)); //re-bind if address is already in use

  memset(&servaddr, 0, sizeof(servaddr));
  memset(&cliaddr, 0, sizeof(cliaddr));
  
  // Filling server information
  servaddr.sin_family    = AF_INET; // IPv4
  servaddr.sin_port = htons(PORT);
  servaddr.sin_addr.s_addr = inet_addr("130.233.123.190");
  if (franka3){
      servaddr.sin_addr.s_addr = inet_addr("130.233.123.182");
  }
  
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

  //Stopper that waits for single msg
  while(true){
      n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
      0, ( struct sockaddr *) &cliaddr,
      &len);

      std::ofstream output(outfile); //ADDED for file-based stopper
      std::cout << "stop message received" << std::endl;
  }


  //receive first message, no timeout for this message
  // n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
  // 0, ( struct sockaddr *) &cliaddr,
  // &len);
  // std::cout << "first message received" << std::endl;

  // //set timeout for subsequent messages
  // struct timeval tv;
  // tv.tv_sec = 0; //timeout_in_seconds
  // tv.tv_usec = 3000; //timeout in microseconds, now 3ms
  // setsockopt(sockfd, SOL_SOCKET, SO_RCVTIMEO, (const char*)&tv, sizeof tv);

  // bool reset = 1;

  // while(1){
  //   n = recvfrom(sockfd, (char *)buffer, MAXLINE,  
  //   0, ( struct sockaddr *) &cliaddr,
  //   &len);

  //   if(buffer[0]=='1'){
  //     std::cout << "reset received" << std::endl;
  //     //remove(outfile.c_str());
  //     reset = 1;
  //   }

  //   else if((n <= 0) && (reset == 1)){
  //     std::ofstream output(outfile); //ADDED for file-based stopper
  //     std::cout << "stopped" << std::endl;
  //     reset = 0;
  //   }

  // }

  return 0;
}