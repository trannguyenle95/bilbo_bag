/* 
Copyright © 2012 NaturalPoint Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. */


#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <fstream> //NOTE OWN ADDITION
#include <ostream> //NOTE OWN ADDITION
#include <iostream> //NOTE OWN ADDITION

#include <vector>
#include <NatNetTypes.h>
#include <NatNetCAPI.h>
#include <NatNetClient.h>

#ifdef _WIN32
#   include <conio.h>
#else
#   include <unistd.h>
#   include <termios.h>
#endif

void NATNET_CALLCONV ServerDiscoveredCallback( const sNatNetDiscoveredServer* pDiscoveredServer, void* pUserContext );
void NATNET_CALLCONV DataHandler(sFrameOfMocapData* data, void* pUserData);    // receives data from the server
void NATNET_CALLCONV MessageHandler(Verbosity msgType, const char* msg);      // receives NatNet error messages
int ConnectClient();

static const ConnectionType kDefaultConnectionType = ConnectionType_Multicast;

NatNetClient* g_pClient = NULL;

std::vector< sNatNetDiscoveredServer > g_discoveredServers;
sNatNetClientConnectParams g_connectParams;
char g_discoveredMulticastGroupAddr[kNatNetIpv4AddrStrLenMax] = NATNET_DEFAULT_MULTICAST_ADDRESS;
int g_analogSamplesPerMocapFrame = 0;
sServerDescription g_serverDescription;


int STOP = 0;

int main( int argc, char* argv[] )
{
    sleep(10); //sleep 10s to give time for setting bag position

    // Install logging callback
    NatNet_SetLogCallback( MessageHandler );

    // create NatNet client
    g_pClient = new NatNetClient();

    // set the frame callback handler
    g_pClient->SetFrameReceivedCallback( DataHandler, g_pClient );	// this function will receive data from the server

    // If no arguments were specified on the command line,
    // attempt to discover servers on the local network.

    // An example of synchronous server discovery.
    // const unsigned int kDiscoveryWaitTimeMillisec = 5 * 1000; // Wait 5 seconds for responses.
    // const int kMaxDescriptions = 10; // Get info for, at most, the first 10 servers to respond.
    // sNatNetDiscoveredServer servers[kMaxDescriptions];
    // int actualNumDescriptions = kMaxDescriptions;
    // NatNet_BroadcastServerDiscovery( servers, &actualNumDescriptions );

    // if ( actualNumDescriptions < kMaxDescriptions )
    // {
    //     // If this happens, more servers responded than the array was able to store.
    // }


    // Do asynchronous server discovery.
    printf( "Looking for servers on the local network.\n" );
    printf( "Press the number key that corresponds to any discovered server to connect to that server.\n" );
    printf( "Press Q at any time to quit.\n\n" );

    NatNetDiscoveryHandle discovery;
    NatNet_CreateAsyncServerDiscovery( &discovery, ServerDiscoveredCallback );

    printf( "before nput \n" );
    
    sleep(1); //Add sleep here so it has time to discover server before trying to connect (?)

    const size_t serverIndex = 0;
    const sNatNetDiscoveredServer& discoveredServer = g_discoveredServers[serverIndex];

    if ( discoveredServer.serverDescription.bConnectionInfoValid )
    {
        // Build the connection parameters.
        g_connectParams.connectionType = discoveredServer.serverDescription.ConnectionMulticast ? ConnectionType_Multicast : ConnectionType_Unicast;
        g_connectParams.serverCommandPort = discoveredServer.serverCommandPort;
        g_connectParams.serverDataPort = discoveredServer.serverDescription.ConnectionDataPort;
        g_connectParams.serverAddress = discoveredServer.serverAddress;
        g_connectParams.localAddress = discoveredServer.localAddress;
        g_connectParams.multicastAddress = g_discoveredMulticastGroupAddr;
    }
    else
    {
        // We're missing some info because it's a legacy server.
        // Guess the defaults and make a best effort attempt to connect.
        g_connectParams.connectionType = kDefaultConnectionType;
        g_connectParams.serverCommandPort = discoveredServer.serverCommandPort;
        g_connectParams.serverDataPort = 0;
        g_connectParams.serverAddress = discoveredServer.serverAddress;
        g_connectParams.localAddress = discoveredServer.localAddress;
        g_connectParams.multicastAddress = NULL;
    }

    NatNet_FreeAsyncServerDiscovery( discovery );

    int iResult;

    // Connect to Motive
    iResult = ConnectClient();
    if (iResult != ErrorCode_OK)
    {
        printf("Error initializing client. See log for details. Exiting.\n");
        return 1;
    }
    else
    {
        printf("Client initialized and ready.\n");
    }
    
	// Ready to receive marker stream!
	printf("\nClient is connected to server and listening for data...\n");
	bool bExit = false;


    int i = 1;
    while(STOP == 0){
        i++;
    }

	// Done - clean up.
	if (g_pClient)
	{
		g_pClient->Disconnect();
		delete g_pClient;
		g_pClient = NULL;
	}

	return ErrorCode_OK;
}


void NATNET_CALLCONV ServerDiscoveredCallback( const sNatNetDiscoveredServer* pDiscoveredServer, void* pUserContext )
{
    char serverHotkey = '.';
    if ( g_discoveredServers.size() < 9 )
    {
        serverHotkey = static_cast<char>('1' + g_discoveredServers.size());
    }

    printf( "[%c] %s %d.%d at %s ",
        serverHotkey,
        pDiscoveredServer->serverDescription.szHostApp,
        pDiscoveredServer->serverDescription.HostAppVersion[0],
        pDiscoveredServer->serverDescription.HostAppVersion[1],
        pDiscoveredServer->serverAddress );

    if ( pDiscoveredServer->serverDescription.bConnectionInfoValid )
    {
        printf( "(%s)\n", pDiscoveredServer->serverDescription.ConnectionMulticast ? "multicast" : "unicast" );
    }
    else
    {
        printf( "(WARNING: Legacy server, could not autodetect settings. Auto-connect may not work reliably.)\n" );
    }

    g_discoveredServers.push_back( *pDiscoveredServer );
}

// Establish a NatNet Client connection
int ConnectClient()
{
    // Release previous server
    g_pClient->Disconnect();

    // Init Client and connect to NatNet server
    int retCode = g_pClient->Connect( g_connectParams );
    if (retCode != ErrorCode_OK)
    {
        printf("Unable to connect to server.  Error code: %d. Exiting.\n", retCode);
        return ErrorCode_Internal;
    }
    else
    {
        // connection succeeded

        void* pResult;
        int nBytes = 0;
        ErrorCode ret = ErrorCode_OK;

        // print server info
        memset( &g_serverDescription, 0, sizeof( g_serverDescription ) );
        ret = g_pClient->GetServerDescription( &g_serverDescription );
        if ( ret != ErrorCode_OK || ! g_serverDescription.HostPresent )
        {
            printf("Unable to connect to server. Host not present. Exiting.\n");
            return 1;
        }
        printf("\n[SampleClient] Server application info:\n");
        printf("Application: %s (ver. %d.%d.%d.%d)\n", g_serverDescription.szHostApp, g_serverDescription.HostAppVersion[0],
            g_serverDescription.HostAppVersion[1], g_serverDescription.HostAppVersion[2], g_serverDescription.HostAppVersion[3]);
        printf("NatNet Version: %d.%d.%d.%d\n", g_serverDescription.NatNetVersion[0], g_serverDescription.NatNetVersion[1],
            g_serverDescription.NatNetVersion[2], g_serverDescription.NatNetVersion[3]);
        printf("Client IP:%s\n", g_connectParams.localAddress );
        printf("Server IP:%s\n", g_connectParams.serverAddress );
        printf("Server Name:%s\n", g_serverDescription.szHostComputerName);

        // get mocap frame rate
        ret = g_pClient->SendMessageAndWait("FrameRate", &pResult, &nBytes);
        if (ret == ErrorCode_OK)
        {
            float fRate = *((float*)pResult);
            printf("Mocap Framerate : %3.2f\n", fRate);
        }
        else
            printf("Error getting frame rate.\n");

        // get # of analog samples per mocap frame of data
        ret = g_pClient->SendMessageAndWait("AnalogSamplesPerMocapFrame", &pResult, &nBytes);
        if (ret == ErrorCode_OK)
        {
            g_analogSamplesPerMocapFrame = *((int*)pResult);
            printf("Analog Samples Per Mocap Frame : %d\n", g_analogSamplesPerMocapFrame);
        }
        else
            printf("Error getting Analog frame rate.\n");
    }

    return ErrorCode_OK;
}

// DataHandler receives data from the server
// This function is called by NatNet when a frame of mocap data is available
void NATNET_CALLCONV DataHandler(sFrameOfMocapData* data, void* pUserData)
{
    NatNetClient* pClient = (NatNetClient*) pUserData;

    int i=0;

    printf("FrameID : %d\n", data->iFrame);
    printf("Timestamp : %3.2lf\n", data->fTimestamp);

	// labeled markers - this includes all markers (Active, Passive, and 'unlabeled' (markers with no asset but a PointCloud ID)
    bool bOccluded;     // marker was not visible (occluded) in this frame
    bool bPCSolved;     // reported position provided by point cloud solve
    bool bModelSolved;  // reported position provided by model solve
    bool bHasModel;     // marker has an associated asset in the data stream
    bool bUnlabeled;    // marker is 'unlabeled', but has a point cloud ID that matches Motive PointCloud ID (In Motive 3D View)
	bool bActiveMarker; // marker is an actively labeled LED marker

	printf("Markers [Count=%d]\n", data->nLabeledMarkers);

    //OWN ADDITION to write most recent (overwrite each time) values to csv:
    std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
    std::ofstream out_file(HOME + "/catkin_ws/src/Data/tracked_markers.csv", std::ios::trunc); //use ios::trunc to OVERWRITE
    if (!out_file)
    {
        std::cout << "Error in creating file!!!" << std::endl; //does file creation fail?
    }

	for(i=0; i < data->nLabeledMarkers; i++)
	{
        bOccluded = ((data->LabeledMarkers[i].params & 0x01)!=0);
        bPCSolved = ((data->LabeledMarkers[i].params & 0x02)!=0);
        bModelSolved = ((data->LabeledMarkers[i].params & 0x04) != 0);
        bHasModel = ((data->LabeledMarkers[i].params & 0x08) != 0);
        bUnlabeled = ((data->LabeledMarkers[i].params & 0x10) != 0);
		bActiveMarker = ((data->LabeledMarkers[i].params & 0x20) != 0);

        sMarker marker = data->LabeledMarkers[i];

        // Marker ID Scheme:
        // Active Markers:
        //   ID = ActiveID, correlates to RB ActiveLabels list
        // Passive Markers: 
        //   If Asset with Legacy Labels
        //      AssetID 	(Hi Word)
        //      MemberID	(Lo Word)
        //   Else
        //      PointCloud ID
        int modelID, markerID;
        NatNet_DecodeID( marker.ID, &modelID, &markerID );
		
        char szMarkerType[512];
        if (bActiveMarker)
            strcpy(szMarkerType, "Active");
        else if(bUnlabeled)
            strcpy(szMarkerType, "Unlabeled");
        else
            strcpy(szMarkerType, "Labeled");

        printf("%s Marker [ModelID=%d, MarkerID=%d] [size=%3.2f] [pos=%3.2f,%3.2f,%3.2f]\n",
            szMarkerType, modelID, markerID, marker.size, marker.x, marker.y, marker.z);

        //OWN ADDITION to write most recent (overwrite each time) values to csv:
        out_file << marker.x << "," << marker.y << "," << marker.z << std::endl;

	}

    STOP = 1; //OWN ADDITION TO MAKE IT STOP AFTER ONE READ
    
}


// MessageHandler receives NatNet error/debug messages
void NATNET_CALLCONV MessageHandler( Verbosity msgType, const char* msg )
{
    // Optional: Filter out debug messages
    if ( msgType < Verbosity_Info )
    {
        return;
    }

    printf( "\n[NatNetLib]" );

    switch ( msgType )
    {
        case Verbosity_Debug:
            printf( " [DEBUG]" );
            break;
        case Verbosity_Info:
            printf( "  [INFO]" );
            break;
        case Verbosity_Warning:
            printf( "  [WARN]" );
            break;
        case Verbosity_Error:
            printf( " [ERROR]" );
            break;
        default:
            printf( " [?????]" );
            break;
    }

    printf( ": %s\n", msg );
}
