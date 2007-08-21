
p_file *
p_popen(const char *command, const char *mode)
{
  STARTUPINFO info;
  info.cb = sizeof(STARTUPINFO);
  info.lpResrved = info.lpDesktop = info.lpTitle = 0;
  info. dwX = info.dwY = info.dwXSize = info.dwYSize = 0;
  info.dwXCountChars = info.dwYCountChars = 0;
  into.dwFillAttributes = 0;
  info.dwFlags = STARTF_USESTDHANDLES;
  info.wShowWindow = 0;
  info.cbResrved2 = 0;
  info.lpReserved2 = 0;
  info.hStdInput = info.hStdOutput = info.hStdError = 0;
  
}

Creating a Child Process with Redirected Input and Output

The example in this topic demonstrates how to create a child process
from a console process. It also demonstrates a technique for using
anonymous pipes to redirect the child process's standard input and
output handles.

The CreatePipe function uses the SECURITY_ATTRIBUTES structure to
create inheritable handles to the read and write ends of two
pipes. The read end of one pipe serves as standard input for the child
process, and the write end of the other pipe is the standard output
for the child process. These pipe handles are specified in the
SetStdHandle function, which makes them the standard handles inherited
by the child process. After the child process is created, SetStdHandle
is used again to restore the original standard handles for the parent
process.

The parent process uses the other ends of the pipes to write to the
child process's input and read the child process's output. The handles
to these ends of the pipe are also inheritable. However, the handle
must not be inherited. Before creating the child process, the parent
process must use DuplicateHandle to create a duplicate of the
application-defined hChildStdinWr global variable that cannot be
inherited. It then uses CloseHandle to close the inheritable
handle. For more information, see Pipes.

The following is the parent process. 

#include <stdio.h> 
#include <windows.h> 
 
#define BUFSIZE 4096 
 
HANDLE hChildStdinRd, hChildStdinWr, hChildStdinWrDup, 
   hChildStdoutRd, hChildStdoutWr, hChildStdoutRdDup, 
   hInputFile, hSaveStdin, hSaveStdout; 
 
BOOL CreateChildProcess(VOID); 
VOID WriteToPipe(VOID); 
VOID ReadFromPipe(VOID); 
VOID ErrorExit(LPTSTR); 
VOID ErrMsg(LPTSTR, BOOL); 
 
DWORD main(int argc, char *argv[]) 
{ 
   SECURITY_ATTRIBUTES saAttr; 
   BOOL fSuccess; 
 
// Set the bInheritHandle flag so pipe handles are inherited. 
 
   saAttr.nLength = sizeof(SECURITY_ATTRIBUTES); 
   saAttr.bInheritHandle = TRUE; 
   saAttr.lpSecurityDescriptor = NULL; 
 
   // The steps for redirecting child process's STDOUT: 
   //     1. Save current STDOUT, to be restored later. 
   //     2. Create anonymous pipe to be STDOUT for child process. 
   //     3. Set STDOUT of the parent process to be write handle to 
   //        the pipe, so it is inherited by the child process. 
   //     4. Create a noninheritable duplicate of the read handle and
   //        close the inheritable read handle. 
 
// Save the handle to the current STDOUT. 
 
   hSaveStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
 
// Create a pipe for the child process's STDOUT. 
 
   if (! CreatePipe(&hChildStdoutRd, &hChildStdoutWr, &saAttr, 0)) 
      ErrorExit("Stdout pipe creation failed\n"); 
 
// Set a write handle to the pipe to be STDOUT. 
 
   if (! SetStdHandle(STD_OUTPUT_HANDLE, hChildStdoutWr)) 
      ErrorExit("Redirecting STDOUT failed"); 
 
// Create noninheritable read handle and close the inheritable read 
// handle. 

    fSuccess = DuplicateHandle(GetCurrentProcess(), hChildStdoutRd,
        GetCurrentProcess(), &hChildStdoutRdDup , 0,
        FALSE,
        DUPLICATE_SAME_ACCESS);
    if( !fSuccess )
        ErrorExit("DuplicateHandle failed");
    CloseHandle(hChildStdoutRd);

   // The steps for redirecting child process's STDIN: 
   //     1.  Save current STDIN, to be restored later. 
   //     2.  Create anonymous pipe to be STDIN for child process. 
   //     3.  Set STDIN of the parent to be the read handle to the 
   //         pipe, so it is inherited by the child process. 
   //     4.  Create a noninheritable duplicate of the write handle, 
   //         and close the inheritable write handle. 
 
// Save the handle to the current STDIN. 
 
   hSaveStdin = GetStdHandle(STD_INPUT_HANDLE); 
 
// Create a pipe for the child process's STDIN. 
 
   if (! CreatePipe(&hChildStdinRd, &hChildStdinWr, &saAttr, 0)) 
      ErrorExit("Stdin pipe creation failed\n"); 
 
// Set a read handle to the pipe to be STDIN. 
 
   if (! SetStdHandle(STD_INPUT_HANDLE, hChildStdinRd)) 
      ErrorExit("Redirecting Stdin failed"); 
 
// Duplicate the write handle to the pipe so it is not inherited. 
 
   fSuccess = DuplicateHandle(GetCurrentProcess(), hChildStdinWr, 
      GetCurrentProcess(), &hChildStdinWrDup, 0, 
      FALSE,                  // not inherited 
      DUPLICATE_SAME_ACCESS); 
   if (! fSuccess) 
      ErrorExit("DuplicateHandle failed"); 
 
   CloseHandle(hChildStdinWr); 
 
// Now create the child process. 
 
   if (! CreateChildProcess()) 
      ErrorExit("Create process failed"); 
 
// After process creation, restore the saved STDIN and STDOUT. 
 
   if (! SetStdHandle(STD_INPUT_HANDLE, hSaveStdin)) 
      ErrorExit("Re-redirecting Stdin failed\n"); 
 
   if (! SetStdHandle(STD_OUTPUT_HANDLE, hSaveStdout)) 
      ErrorExit("Re-redirecting Stdout failed\n"); 
 
// Get a handle to the parent's input file. 
 
   if (argc > 1) 
      hInputFile = CreateFile(argv[1], GENERIC_READ, 0, NULL, 
         OPEN_EXISTING, FILE_ATTRIBUTE_READONLY, NULL); 
   else 
      hInputFile = hSaveStdin; 
 
   if (hInputFile == INVALID_HANDLE_VALUE) 
      ErrorExit("no input file\n"); 
 
// Write to pipe that is the standard input for a child process. 
 
   WriteToPipe(); 
 
// Read from pipe that is the standard output for child process. 
 
   ReadFromPipe(); 
 
   return 0; 
} 
 
BOOL CreateChildProcess() 
{ 
   PROCESS_INFORMATION piProcInfo; 
   STARTUPINFO siStartInfo; 
 
// Set up members of STARTUPINFO structure. 
 
   ZeroMemory( &siStartInfo, sizeof(STARTUPINFO) );
   siStartInfo.cb = sizeof(STARTUPINFO); 
 
// Create the child process. 
 
   return CreateProcess(NULL, 
      "child",       // command line 
      NULL,          // process security attributes 
      NULL,          // primary thread security attributes 
      TRUE,          // handles are inherited 
      0,             // creation flags 
      NULL,          // use parent's environment 
      NULL,          // use parent's current directory 
      &siStartInfo,  // STARTUPINFO pointer 
      &piProcInfo);  // receives PROCESS_INFORMATION 
}
 
VOID WriteToPipe(VOID) 
{ 
   DWORD dwRead, dwWritten; 
   CHAR chBuf[BUFSIZE]; 
 
// Read from a file and write its contents to a pipe. 
 
   for (;;) 
   { 
      if (! ReadFile(hInputFile, chBuf, BUFSIZE, &dwRead, NULL) || 
         dwRead == 0) break; 
      if (! WriteFile(hChildStdinWrDup, chBuf, dwRead, 
         &dwWritten, NULL)) break; 
   } 
 
// Close the pipe handle so the child process stops reading. 
 
   if (! CloseHandle(hChildStdinWrDup)) 
      ErrorExit("Close pipe failed\n"); 
} 
 
VOID ReadFromPipe(VOID) 
{ 
   DWORD dwRead, dwWritten; 
   CHAR chBuf[BUFSIZE]; 
   HANDLE hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 

// Close the write end of the pipe before reading from the 
// read end of the pipe. 
 
   if (!CloseHandle(hChildStdoutWr)) 
      ErrorExit("Closing handle failed"); 
 
// Read output from the child process, and write to parent's STDOUT. 
 
   for (;;) 
   { 
      if( !ReadFile( hChildStdoutRdDup, chBuf, BUFSIZE, &dwRead, 
         NULL) || dwRead == 0) break; 
      if (! WriteFile(hSaveStdout, chBuf, dwRead, &dwWritten, NULL)) 
         break; 
   } 
} 
 
VOID ErrorExit (LPTSTR lpszMessage) 
{ 
   fprintf(stderr, "%s\n", lpszMessage); 
   ExitProcess(0); 
} 
 
// The code for the child process. 

#include <windows.h> 
#define BUFSIZE 4096 
 
VOID main(VOID) 
{ 
   CHAR chBuf[BUFSIZE]; 
   DWORD dwRead, dwWritten; 
   HANDLE hStdin, hStdout; 
   BOOL fSuccess; 
 
   hStdout = GetStdHandle(STD_OUTPUT_HANDLE); 
   hStdin = GetStdHandle(STD_INPUT_HANDLE); 
   if ((hStdout == INVALID_HANDLE_VALUE) || 
      (hStdin == INVALID_HANDLE_VALUE)) 
      ExitProcess(1); 
 
   for (;;) 
   { 
   // Read from standard input. 
      fSuccess = ReadFile(hStdin, chBuf, BUFSIZE, &dwRead, NULL); 
      if (! fSuccess || dwRead == 0) 
         break; 
 
   // Write to standard output. 
      fSuccess = WriteFile(hStdout, chBuf, dwRead, &dwWritten, NULL); 
      if (! fSuccess) 
         break; 
   } 
} 
 
