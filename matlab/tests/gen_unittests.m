% Generates mat files for loadmat unit tests
% This is the version for matlab 5 and higher
% Uses save_test.m function

% work out matlab version and file suffix for test files
global FILEPREFIX FILESUFFIX
FILEPREFIX = [fullfile(pwd, 'data') filesep];
temp = ver('MATLAB');
mlv = temp.Version;
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


% struct
save_test('teststruct', ...
	  struct('stringfield','Rats live on no evil star.',...
		 'doublefield',[sqrt(2) exp(1) pi],...
		 'complexfield',(1+1j)*[sqrt(2) exp(1) pi]));

% cell
save_test('testcell', ...
	  {['This cell contains this string and 3 arrays of increasing' ...
	    ' length'], 1., 1.:2., 1.:3.});

% Empty cells in two cell matrices
save_test('testemptycell', {1, 2, [], [], 3});

% 3D matrix
save_test('test3dmatrix', reshape(1:24,[2 3 4]))

% nested cell array
save_test('testcellnest', {1, {2, 3, {4, 5}}});

% nested struct
save_test('teststructnest', struct('one', 1, 'two', ...
				   struct('three', 'number 3')));

% array of struct
save_test('teststructarr', [struct('one', 1, 'two', 2) ...
		    struct('one', 'number 1', 'two', 'number 2')]);

% matlab object
save_test('testobject', inline('x'))

% array of matlab objects
%save_test('testobjarr', [inline('x') inline('y')])

% unicode test
if str2num(mlv) > 7  % function added 7.0.1
  fid = fopen([FILEPREFIX 'japanese_utf8.txt']);
  from_japan = fread(fid, 'uint8')';
  fclose(fid);
  save_test('testunicode', native2unicode(from_japan, 'utf-8'));
end
  