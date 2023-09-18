/* Definitions for use with Linux AF_PACKET sockets.
   Copyright (C) 1998, 1999 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, write to the Free
   Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
   02111-1307 USA.  */

#ifndef __NETPACKET_PACKET_H
#define __NETPACKET_PACKET_H	1

struct sockaddr_ll
  {
    unsigned short int sll_family;
    unsigned short int sll_protocol;
    int sll_ifindex;
    unsigned short int sll_hatype;
    unsigned char sll_pkttype;
    unsigned char sll_halen;
    unsigned char sll_addr[8];
  };

/* Packet types.  */

#define PACKET_HOST		0		/* To us.  */
#define PACKET_BROADCAST	1		/* To all.  */
#define PACKET_MULTICAST	2		/* To group.  */
#define PACKET_OTHERHOST	3		/* To someone else.  */
#define PACKET_OUTGOING		4		/* Originated by us . */
#define PACKET_LOOPBACK		5
#define PACKET_FASTROUTE	6

/* Packet socket options.  */

#define PACKET_ADD_MEMBERSHIP		1
#define PACKET_DROP_MEMBERSHIP		2
#define	PACKET_RECV_OUTPUT		3
#define	PACKET_RX_RING			5
#define	PACKET_STATISTICS		6

struct packet_mreq
  {
    int mr_ifindex;
    unsigned short int mr_type;
    unsigned short int mr_alen;
    unsigned char mr_address[8];
  };

#define PACKET_MR_MULTICAST	0
#define PACKET_MR_PROMISC	1
#define PACKET_MR_ALLMULTI	2

#endif	/* netpacket/packet.h */
