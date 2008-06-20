% Generates mat files for loadmat unit tests
% Uses save_test.m function
% This is the version for matlab 4

% work out matlab version and file suffix for test files
global FILEPREFIX FILESUFFIX
sepchar = '/';
if strcmp(computer, 'PCWIN'), sepchar = '\'; end
FILEPREFIX = [pwd sepchar 'data' sepchar];
mlv = version;
FILESUFFIX = ['_' mlv '_' computer '.mat'];

% basic double array
save_test('testdouble', 0:pi/4:2*pi);

% string
save_test('teststring', '"Do nine men interpret?" "Nine men," I nod.')

% complex
theta = 0:pi/4:2*pi;
save_test('testcomplex', cos(theta) + 1j*sin(theta));

% asymmetric array to check indexing
a = zeros(3, 5);
a(:,1) = [1:3]';
a(1,:) = 1:5;

% 2D matrix
save_test('testmatrix', a);

% minus number - tests signed int 
save_test('testminus', -1);

% single character
save_test('testonechar', 'r');

% string array
save_test('teststringarray', ['one  '; 'two  '; 'three']);

% sparse array
save_test('testsparse', sparse(a));

% sparse complex array
b = sparse(a);
b(1,1) = b(1,1) + j;
save_test('testsparsecomplex', b);

% Two variables in same file
save([FILEPREFIX 'testmulti' FILESUFFIX], 'a', 'theta')

