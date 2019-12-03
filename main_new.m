%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ���C�����[�`��
% (1)�e�X�g�s��쐬
% (2)FORTRAN�ł�dqds�̎��s
% (3)MATLAB�ł�dqds�̎��s
% (4)���{�����x���Z�ɂ��^�l�̑�p���̌v�Z
% (5)�덷�]��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�@���̃t�H���_���ɕK�v�ȃt�@�C��
% main.m 
% dqds_2_1_6_0_0_test.m 
% test_dpteqr.f90 
% dpteqr_O3.out 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORTRAN�ŃR���p�C�����āCtest_dpteqr.f90����dpteqr_O3.out�𐶐�
% �ȉ��̃R�}���h���^�[�~�i���Ŏ��s�F
% gfortran -Wall -O3 test_dpteqr.f90 -o dpteqr /usr/local/lib/liblapack.a /usr/local/lib/librefblas.a 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)�e�X�g�s��쐬
% �y���͍s��zA2
%  matrix14�̃e�X�g�s��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('(1)�e�X�g�s��쐬');
rand('seed',1);
% �s��̃T�C�Y���w�肷��
prompt = '�s��̃T�C�Y���w�肵�Ă������� ';
x = input(prompt)
n=int32(x);
% �e�X�g�s��̔ԍ����w�肷��
prompt = '�e�X�g�s��̔ԍ���1~20�Ŏw�肵�Ă��������@';
y = input(prompt)

% test_matirix.m�̊֐����Ăяo��
[A,a,b]=test_matrix(n,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)FORTRAN�ł�dqds�̎��s
% �y�o�͌��ʁz�ŗL�l�Feig_dpteqr�C�v�Z���ԁFtime_dpteqr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(2)FORTRAN�ł�dqds�̎��s');
% �s��A2���o�C�i���`���̃t�@�C��(�t�@�C�����Fin.bin)�Ƃ��ĕۑ�
dmat_bin_save(A,'in.bin');
dvec_bin_save(b,'in_diag1.bin');
dvec_bin_save(a,'in_diag2.bin');
% type='DST';
% fileID = fopen('in.bin','w');
% fwrite(fileID,type,'char','ieee-le');
% fwrite(fileID,n,'int32','ieee-le');
% fwrite(fileID,A,'double','ieee-le');
% fclose(fileID);


% FORTRAN�ł�dqds�̎��s�D�O���R�}���hdpteqr�����s�iLAPACK���C�u������dpteqr���[�`���j
command='./dpteqr';
disp(sprintf('%s',command));
[status,cmdout]=system(command);
disp(sprintf('statud=%d',status));
disp(sprintf('%s',cmdout));

% �o�C�i���`���̏o�͌��ʂ̃t�@�C��(�t�B�C�����Fout_dpteqr.bin�j��ǂݍ���
fileID=fopen('out_dpteqr.bin','r');status = fseek(fileID,0,'bof');
type_dpteqr=fread(fileID,4,'*char','ieee-le');type_dpteqr=type_dpteqr';
n=fread(fileID,1,'*int32','ieee-le');A_dpteqr= fread(fileID,[n n],'double','ieee-le');
%n1�ɕύX���遪
routine_dpteqr=fread(fileID,15,'*char','ieee-le');routine_dpteqr=routine_dpteqr';
% �v�Z���Ԃ̓ǂݍ���
time_dpteqr=fread(fileID,1,'double','ieee-le');
% �o�͌��ʂ̌ŗL�l��ǂݍ��݃\�[�g
eig_dpteqr_r=fread(fileID,[n 1],'double','ieee-le');
eig_dpteqr=sort(eig_dpteqr_r);
fclose(fileID);
disp(sprintf('Computing time=%g [s]',time_dpteqr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% C dqds�̌��ʂ�ǂݍ���       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(2.5)C�ł�dqds�̎��s');

command='./dqds_double';
disp(sprintf('%s',command));
[status,cmdout]=system(command);    
dqds_c_in=dvec_bin_load('out_c_dqds.bin');
% �����񐔂̓ǂݍ���
dqds_c_times_in=ivec_load('repeat_times.txt');
%�\�[�g
[dqds_c,I]=sort(dqds_c_in);
dqds_c_times=dqds_c_times_in(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)MATLAB�ł�dqds�̎��s
% �y�o�͌��ʁz�ŗL�l�Flambda3�C�v�Z���ԁFtime_matlab_dqds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(3)MATLAB�ł�dqds�̎��s');
tic;t=toc;
% dqds���s
tic;
debug=1;
[lambda3,nn3]=dqds_2_1_6_0_0_test(n,a,b,debug);
time_matlab_dqds=toc;
% �v�Z����
disp(sprintf('Computing time=%g [s]',time_matlab_dqds));
% dqds�ł̌ŗL�l���\�[�g
[lambda3,I3]=sort(lambda3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4)���{�����x���Z�ɂ��^�l�̑�p���̌v�Z
% �y�o�͌��ʁz�ŗL�l�Flambda�C�v�Z���ԁFtime_multi_eig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(4)���{�����x���Z�ɂ��^�l�̑�p���̌v�Z');
% �g�p���� prec [bits]
set_default_prec(512);
%�@���{�����x���Z�ŌŗL�l�v�Z
tic;
lambda=double(eig(multi(A)));
time_multi_eig=toc;
disp(sprintf('Computing time=%g [s]',time_multi_eig));
% �{���x�ɃL���X�g���C�\�[�g���Ă���C�܂����{���ɃL���X�g
lambda=sort(double(lambda));
lambda=double(lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5)�덷�]��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(5)�덷�]��');
erzero=ones(n,1).*2^(-53);


% MATLAB�ł�dqds�̌덷
% MATLAB�̌덷
disp(sprintf('=== dqds (MATLAB)'));
er3=double((abs((lambda3-lambda)./lambda)));
er3zero=max(er3,erzero);
er3_max=max(er3zero);
er3_ave=10^(sum(log10(er3zero))/double(n));
disp(sprintf('Max: %.1e',er3_max));
disp(sprintf('Ave: %.1e',er3_ave));

% % FORTRAN�ł�dqds�̌덷
disp(sprintf('=== dpteqr (FORTRAN)'));
lambda4=eig_dpteqr;
er4=double(abs((lambda4-lambda)./lambda));
er4zero=max(er4,erzero);
er4_max=max(er4zero);
er4_ave=10^(sum(log10(er4zero))/double(n));
disp(sprintf('Max: %.1e',er4_max));
disp(sprintf('Ave: %.1e',er4_ave));

% C�ł�dqds�̌덷
disp(sprintf('=== dqds ( C )'));
lambda5=dqds_c;
er5=double(abs((lambda5-lambda)./lambda));
er5zero=max(er5,erzero);
er5_max=max(er5zero);
er5_ave=10^(sum(log10(er5zero))/double(n));
disp(sprintf('Max: %.1e',er5_max));
disp(sprintf('Ave: %.1e',er5_ave));

disp('(5*)�덷�����������ƈႤ�������i�[');

% �قȂ�덷�̔z����쐬
% MATLAB��FORTRAN
% �덷�������������P�ňႤ�������O�ŏo��
eq = abs(er3zero-er4zero)== 0;
% �덷���Ⴄ�����݂̂��o��
er3zero_not_eq = er3zero(~eq);
er4zero_not_eq = er4zero(~eq);
n1 = numel(er3zero_not_eq);
disp(sprintf('MATLAB��FORTRAN�̌덷���قȂ鐔: %d ��',n1));
% �덷�����������̔z��̂ݏo��
er_zero_eq = er3zero(eq);
disp(sprintf('MATLAB��FORTRAN�̌덷��������: %d ��',numel(er_zero_eq)));

% C��FORTRAN
% �덷�������������P�ňႤ�������O�ŏo��
eq_1 = abs(er5zero-er4zero)== 0;
% �덷���Ⴄ�����݂̂��o��
er5zero_not_eq = er5zero(~eq_1);
er4zero_not_eq_1 = er4zero(~eq_1);
n2 = numel(er5zero_not_eq);
disp(sprintf('C��FORTRAN�̌덷���قȂ鐔: %d ��',n2));
% �덷�����������̔z��̂ݏo��
er_zero_eq_1 = er5zero(eq_1);
disp(sprintf('C��FORTRAN�̌덷��������: %d ��',numel(er_zero_eq)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5*)���ςƕ��U
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('(5*)���U�ƕ���');

% �قȂ�덷�̔z��̕��ςƕ��U
var3 = var(er3zero_not_eq);
var4 = var(er4zero_not_eq);
disp(sprintf('MATLAB�̕��U: %.1e',var3));
disp(sprintf('FORTRAN�̕��U: %.1e',var4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5**)���_�v�Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -14.5��^�񒆂Ƃ��Ă����菬�����Ȃ�v���X�A�傫���Ȃ�}�C�i�X�Ƃ��Čv�Z
disp('(5**)���_��\��');
for i=1:n2
    p3(i) = -14.5-log10(er3zero_not_eq(i));
    p4(i) = -14.5-log10(er4zero_not_eq(i));
end
disp(sprintf('MATLAB�̓��_: %g�_',sum(p3)));
disp(sprintf('FORTRAN�̓��_: %g�_',sum(p4)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6)�O���t
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name_fig = append('matrix',string(x),'_','test',string(y),'_','fig');
mkdir(folder_name_fig);

% �O���t
figure(1);
semilogy(1:double(n),er3zero,'bo',1:double(n),er4zero,'rx',1:double(n),er5zero,'gs','MarkerSize',13,'LineWidth',3,'MarkerFaceColor','w');
grid on;
axis([1 double(n) 1e-17 1e-13]);
xlabel('Eigenvalue number');
ylabel('Relative errors');
legend('dqds (MATLAB)','dpteqr (FORTRAN)','dqds(C)');
saveas(gcf,'all_er.fig')
saveas(gcf,'all_er.png')
movefile('all_er.fig',(folder_name_fig));
movefile('all_er.png',(folder_name_fig));

figure(2);
histogram(log10(er3zero));title('dqds (MATLAB)');axis([-17 -13 0 35]);grid on;
%%%%% �ΐ��̖��O�ǂ�����
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'MATLAB_er.fig')
saveas(gcf,'MATLAB_er.png')
movefile('MATLAB_er.fig',(folder_name_fig));
movefile('MATLAB_er.png',(folder_name_fig));

figure(3);
histogram(log10(er4zero));title('dpteqr (FORTRAN)');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors');
ylabel('Number');
saveas(gcf,'FORTRAN_er.fig')
saveas(gcf,'FORTRAN_er.png')
movefile('FORTRAN_er.fig',(folder_name_fig));
movefile('FORTRAN_er.png',(folder_name_fig));

figure(4);
histogram(log10(er5zero));title('dpteqr (dqds C)');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors');
saveas(gcf,'C_er.fig')
saveas(gcf,'C_er.png')
movefile('C_er.fig',(folder_name_fig));
movefile('C_er.png',(folder_name_fig));

% �قȂ�덷�݂̂��o�͂����O���t(MATLAB��FORTRAN)
% �^�C�g���@���̃q�X�g�O�����H�@���O���߂�@����loge k?

figure(5);
semilogy(1:double(n1),er3zero_not_eq,'bo',1:double(n1),er4zero_not_eq,'rx','MarkerSize',13,'LineWidth',3,'MarkerFaceColor','w');
grid on;
axis([1 double(n1) 1e-17 1e-13]);
ylabel('Relative errors');
legend('dqds (MATLAB)','dpteqr (FORTRAN)');
saveas(gcf,'not_eq_er.fig')
saveas(gcf,'not_eq_er.png')
movefile('not_eq_er.fig',(folder_name_fig));
movefile('not_eq_er.png',(folder_name_fig));

figure(6);
histogram(log10(er5zero_not_eq));title('dqds (C):�덷���قȂ�');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'not_eq_C_er.fig')
saveas(gcf,'not_eq_C_er.png')
movefile('not_eq_C_er.fig',(folder_name_fig));
movefile('not_eq_C_er.png',(folder_name_fig));

figure(7);
histogram(log10(er3zero_not_eq));title('dqds (MATLAB):�덷���قȂ�');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'not_eq_MATLAB_er.fig')
saveas(gcf,'not_eq_MATLAB_er.png')
movefile('not_eq_MATLAB_er.fig',(folder_name_fig));
movefile('not_eq_MATLAB_er.png',(folder_name_fig));

figure(8);
histogram(log10(er4zero_not_eq));title('dqteqr (FORTRAN):�덷���قȂ�');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'not_eq_FORTRAN_er.fig')
saveas(gcf,'not_eq_FORTRAN_er.png')
movefile('not_eq_FORTRAN_er.fig',(folder_name_fig));
movefile('not_eq_FORTRAN_er.png',(folder_name_fig));

figure(9);
semilogy(dqds_c_times,er5zero,'bo','MarkerSize',13,'LineWidth',3,'MarkerFaceColor','w');
grid on;
axis([0 700 1e-17 1e-13]);
xlabel('������');
ylabel('Relative errors');
legend('dqds( C )');
saveas(gcf,'loop_er.fig')
saveas(gcf,'loop_er.png')
movefile('loop_er.fig',(folder_name_fig));
movefile('loop_er.png',(folder_name_fig));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (7)�t�@�C���ւ̏o��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name = append('matrix',string(x),'_','test',string(y));
mkdir(folder_name);
% ��r���ʁi���s���ԁA�������񐔁A�������́A�ő�l�A���ρA���U�A�������ƈႤ���j
file_name = append('data',string(y),'_','matrix','_',string(x),'.txt');
fid = fopen((file_name),'w');
% �g�p�����e�X�g�s��̏��
fprintf(fid,'(1)���\r\n');
name = append('�e�X�g�s��',string(y));
fprintf(fid,'�e�X�g�s��:%s\n',name);
size = string(x);
fprintf(fid,'�s��T�C�Y:%s\n',size);
fprintf(fid,'=============================================\r\n');
% ���s����
fprintf(fid,'(2)���s����\r\n');
fprintf(fid,'FORTRAN_Computing_time %g [s]\r\n',time_dpteqr);
fprintf(fid,'MATLAB_Computing_time %g [s]\r\n',time_matlab_dqds);
fprintf(fid,'=============================================\r\n');
% ��������
fprintf(fid,'(3)��������\r\n');
fprintf(fid,'MATLAB: %d\r\n',nn3);
fprintf(fid,'=============================================\r\n');
% �덷�]��
fprintf(fid,'(4)�덷�]��\r\n');
fprintf(fid,'=============================================\r\n');
fprintf(fid,'\r\n');
% MATLAB 
fprintf(fid,'=== dqds(MATLAB) \r\n');
fprintf(fid,'Max: %.1e\r\n',er3_max);
fprintf(fid,'Ave: %.1e\r\n',er3_ave);
% FORTRAN
fprintf(fid,'=== dqds(FOTRAN) \r\n');
fprintf(fid,'Max: %.1e\r\n',er4_max);
fprintf(fid,'Ave: %.1e\r\n',er4_ave);
% c
fprintf(fid,'=== dqds(C) \r\n');
fprintf(fid,'Max: %.1e\r\n',er5_max);
fprintf(fid,'Ave: %.1e\r\n',er5_ave);
fprintf(fid,'=============================================\r\n');
% �덷�����������ƈႤ����
fprintf(fid,'(5)�덷���������ƈႤ��\r\n');
fprintf(fid,'�덷���قȂ鐔: %d ��\r\n',n1);
fprintf(fid,'�덷��������: %d ��\r\n',numel(er_zero_eq));
fprintf(fid,'=============================================\r\n');
% ���U
fprintf(fid,'(6)�덷���قȂ镔���̕��U\r\n');
fprintf(fid,'MATLAB�̕��U: %.1e\r\n',var3);
fprintf(fid,'FORTRAN�̕��U: %.1e\r\n',var4);
fprintf(fid,'=============================================\r\n');
% ���_
fprintf(fid,'(7)���_\r\n');
fprintf(fid,'MATLAB�̓��_: %g�_\r\n',sum(p3));
fprintf(fid,'FORTRAN�̓��_: %g�_\r\n',sum(p4));
fprintf(fid,'��̓��_��: %g�_\r\n',abs(sum(p3)-sum(p4)));
fprintf(fid,'=============================================\r\n');

fclose(fid);

% �^�l
file_name1 = 'lambda.txt';
fid = fopen((file_name1),'w');
fprintf(fid,'�ŗL�l(�^�l)\r\n');
fprintf(fid,'%g\n',lambda);
fclose(fid);
% FORTRAN�̌ŗL�l
file_name2 = 'FORTRAN_eig.txt';
fid = fopen((file_name2),'w');
fprintf(fid,'�ŗL�l(FORTRAN)\r\n');
fprintf(fid,'%g\n',eig_dpteqr);
fclose(fid);
% MATLAB�̌ŗL�l
file_name3 = 'MATLAB_eig.txt';
fid = fopen((file_name3),'w');
fprintf(fid,'�ŗL�l(MATLAB)\r\n');
fprintf(fid,'%g\n',lambda3);
fclose(fid);
% C�̌ŗL�l
file_name4 = 'C_eig.txt';
fid = fopen((file_name4),'w');
fprintf(fid,'�ŗL�l(c)\r\n');
fprintf(fid,'%g\n',dqds_c);
fclose(fid);

% �덷MATLAB
file_name5 = 'MATLAB_er.txt';
fid = fopen((file_name5),'w');
fprintf(fid,'�덷(c)\r\n');
fprintf(fid,'%g\n',er3);
fclose(fid);
% �덷FORTRAN
file_name6 = 'FORTRAN_er.txt';
fid = fopen((file_name6),'w');
fprintf(fid,'�덷(FORTRAN)\r\n');
fprintf(fid,'%g\n',er4);
fclose(fid);
% �덷C
file_name7 = 'C_er.txt';
fid = fopen((file_name7),'w');
fprintf(fid,'�덷(c)\r\n');
fprintf(fid,'%g\n',er5);
fclose(fid);

movefile((folder_name_fig),(folder_name));
movefile((file_name),(folder_name));
movefile((file_name1),(folder_name));
movefile((file_name2),(folder_name));
movefile((file_name3),(folder_name));
movefile((file_name4),(folder_name));
movefile((file_name5),(folder_name));
movefile((file_name6),(folder_name));
movefile((file_name7),(folder_name));
