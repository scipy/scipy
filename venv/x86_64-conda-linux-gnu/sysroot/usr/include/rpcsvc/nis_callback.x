%/*
% * Sun RPC is a product of Sun Microsystems, Inc. and is provided for
% * unrestricted use provided that this legend is included on all tape
% * media and as a part of the software program in whole or part.  Users
% * may copy or modify Sun RPC without charge, but are not authorized
% * to license or distribute it to anyone else except as part of a product or
% * program developed by the user or with the express written consent of
% * Sun Microsystems, Inc.
% *
% * SUN RPC IS PROVIDED AS IS WITH NO WARRANTIES OF ANY KIND INCLUDING THE
% * WARRANTIES OF DESIGN, MERCHANTIBILITY AND FITNESS FOR A PARTICULAR
% * PURPOSE, OR ARISING FROM A COURSE OF DEALING, USAGE OR TRADE PRACTICE.
% *
% * Sun RPC is provided with no support and without any obligation on the
% * part of Sun Microsystems, Inc. to assist in its use, correction,
% * modification or enhancement.
% *
% * SUN MICROSYSTEMS, INC. SHALL HAVE NO LIABILITY WITH RESPECT TO THE
% * INFRINGEMENT OF COPYRIGHTS, TRADE SECRETS OR ANY PATENTS BY SUN RPC
% * OR ANY PART THEREOF.
% *
% * In no event will Sun Microsystems, Inc. be liable for any lost revenue
% * or profits or other special, indirect and consequential damages, even if
% * Sun has been advised of the possibility of such damages.
% *
% * Sun Microsystems, Inc.
% * 2550 Garcia Avenue
% * Mountain View, California  94043
% */

/*
 *	nis_callback.x
 *
 *	Copyright (c) 1988-1992 Sun Microsystems Inc
 *	All Rights Reserved.
 */

%#pragma ident	"@(#)nis_callback.x	1.7	94/05/03 SMI"

/*
 * "@(#)zns_cback.x 1.2 90/09/10 Copyr 1990 Sun Micro" 
 *
 * RPCL description of the Callback Service.
 */

#ifdef RPC_HDR
%#include <rpcsvc/nis.h>
#endif
#ifdef RPC_XDR
%#include "nis_clnt.h"
#endif

typedef nis_object	*obj_p;

struct cback_data {
	obj_p		entries<>;	/* List of objects */
};

program CB_PROG {
	version CB_VERS {
		bool	CBPROC_RECEIVE(cback_data) = 1;
		void	CBPROC_FINISH(void) = 2;
		void	CBPROC_ERROR(nis_error) = 3;
	} = 1;
} = 100302;
