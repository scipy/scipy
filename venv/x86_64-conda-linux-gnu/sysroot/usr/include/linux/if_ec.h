/* Definitions for Econet sockets. */

#ifndef __LINUX_IF_EC
#define __LINUX_IF_EC

/* User visible stuff. Glibc provides its own but libc5 folk will use these */

struct ec_addr
{
  unsigned char station;		/* Station number.  */
  unsigned char net;			/* Network number.  */
};

struct sockaddr_ec
{
  unsigned short sec_family;
  unsigned char port;			/* Port number.  */
  unsigned char cb;			/* Control/flag byte.  */
  unsigned char type;			/* Type of message.  */
  struct ec_addr addr;
  unsigned long cookie;
};

#define ECTYPE_PACKET_RECEIVED		0	/* Packet received */
#define ECTYPE_TRANSMIT_STATUS		0x10	/* Transmit completed, 
						   low nibble holds status */

#define ECTYPE_TRANSMIT_OK		1
#define ECTYPE_TRANSMIT_NOT_LISTENING	2
#define ECTYPE_TRANSMIT_NET_ERROR	3
#define ECTYPE_TRANSMIT_NO_CLOCK	4
#define ECTYPE_TRANSMIT_LINE_JAMMED	5
#define ECTYPE_TRANSMIT_NOT_PRESENT	6


#endif
