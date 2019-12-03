%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% メインルーチン
% (1)テスト行列作成
% (2)FORTRAN版のdqdsの実行
% (3)MATLAB版のdqdsの実行
% (4)多倍長精度演算により真値の代用物の計算
% (5)誤差評価
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%　このフォルダ内に必要なファイル
% main.m 
% dqds_2_1_6_0_0_test.m 
% test_dpteqr.f90 
% dpteqr_O3.out 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORTRANでコンパイルして，test_dpteqr.f90からdpteqr_O3.outを生成
% 以下のコマンドをターミナルで実行：
% gfortran -Wall -O3 test_dpteqr.f90 -o dpteqr /usr/local/lib/liblapack.a /usr/local/lib/librefblas.a 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (1)テスト行列作成
% 【入力行列】A2
%  matrix14のテスト行列
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('(1)テスト行列作成');
rand('seed',1);
% 行列のサイズを指定する
prompt = '行列のサイズを指定してください ';
x = input(prompt)
n=int32(x);
% テスト行列の番号を指定する
prompt = 'テスト行列の番号を1~20で指定してください　';
y = input(prompt)

% test_matirix.mの関数を呼び出す
[A,a,b]=test_matrix(n,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (2)FORTRAN版のdqdsの実行
% 【出力結果】固有値：eig_dpteqr，計算時間：time_dpteqr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(2)FORTRAN版のdqdsの実行');
% 行列A2をバイナリ形式のファイル(ファイル名：in.bin)として保存
dmat_bin_save(A,'in.bin');
dvec_bin_save(b,'in_diag1.bin');
dvec_bin_save(a,'in_diag2.bin');
% type='DST';
% fileID = fopen('in.bin','w');
% fwrite(fileID,type,'char','ieee-le');
% fwrite(fileID,n,'int32','ieee-le');
% fwrite(fileID,A,'double','ieee-le');
% fclose(fileID);


% FORTRAN版のdqdsの実行．外部コマンドdpteqrを実行（LAPACKライブラリのdpteqrルーチン）
command='./dpteqr';
disp(sprintf('%s',command));
[status,cmdout]=system(command);
disp(sprintf('statud=%d',status));
disp(sprintf('%s',cmdout));

% バイナリ形式の出力結果のファイル(フィイル名：out_dpteqr.bin）を読み込む
fileID=fopen('out_dpteqr.bin','r');status = fseek(fileID,0,'bof');
type_dpteqr=fread(fileID,4,'*char','ieee-le');type_dpteqr=type_dpteqr';
n=fread(fileID,1,'*int32','ieee-le');A_dpteqr= fread(fileID,[n n],'double','ieee-le');
%n1に変更する↑
routine_dpteqr=fread(fileID,15,'*char','ieee-le');routine_dpteqr=routine_dpteqr';
% 計算時間の読み込み
time_dpteqr=fread(fileID,1,'double','ieee-le');
% 出力結果の固有値を読み込みソート
eig_dpteqr_r=fread(fileID,[n 1],'double','ieee-le');
eig_dpteqr=sort(eig_dpteqr_r);
fclose(fileID);
disp(sprintf('Computing time=%g [s]',time_dpteqr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% C dqdsの結果を読み込む       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(2.5)C版のdqdsの実行');

command='./dqds_double';
disp(sprintf('%s',command));
[status,cmdout]=system(command);    
dqds_c_in=dvec_bin_load('out_c_dqds.bin');
% 反復回数の読み込み
dqds_c_times_in=ivec_load('repeat_times.txt');
%ソート
[dqds_c,I]=sort(dqds_c_in);
dqds_c_times=dqds_c_times_in(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3)MATLAB版のdqdsの実行
% 【出力結果】固有値：lambda3，計算時間：time_matlab_dqds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(3)MATLAB版のdqdsの実行');
tic;t=toc;
% dqds実行
tic;
debug=1;
[lambda3,nn3]=dqds_2_1_6_0_0_test(n,a,b,debug);
time_matlab_dqds=toc;
% 計算時間
disp(sprintf('Computing time=%g [s]',time_matlab_dqds));
% dqdsでの固有値をソート
[lambda3,I3]=sort(lambda3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (4)多倍長精度演算により真値の代用物の計算
% 【出力結果】固有値：lambda，計算時間：time_multi_eig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(4)多倍長精度演算により真値の代用物の計算');
% 使用桁数 prec [bits]
set_default_prec(512);
%　多倍長精度演算で固有値計算
tic;
lambda=double(eig(multi(A)));
time_multi_eig=toc;
disp(sprintf('Computing time=%g [s]',time_multi_eig));
% 倍精度にキャストし，ソートしてから，また多倍長にキャスト
lambda=sort(double(lambda));
lambda=double(lambda);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5)誤差評価
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('(5)誤差評価');
erzero=ones(n,1).*2^(-53);


% MATLAB版のdqdsの誤差
% MATLABの誤差
disp(sprintf('=== dqds (MATLAB)'));
er3=double((abs((lambda3-lambda)./lambda)));
er3zero=max(er3,erzero);
er3_max=max(er3zero);
er3_ave=10^(sum(log10(er3zero))/double(n));
disp(sprintf('Max: %.1e',er3_max));
disp(sprintf('Ave: %.1e',er3_ave));

% % FORTRAN版のdqdsの誤差
disp(sprintf('=== dpteqr (FORTRAN)'));
lambda4=eig_dpteqr;
er4=double(abs((lambda4-lambda)./lambda));
er4zero=max(er4,erzero);
er4_max=max(er4zero);
er4_ave=10^(sum(log10(er4zero))/double(n));
disp(sprintf('Max: %.1e',er4_max));
disp(sprintf('Ave: %.1e',er4_ave));

% C版のdqdsの誤差
disp(sprintf('=== dqds ( C )'));
lambda5=dqds_c;
er5=double(abs((lambda5-lambda)./lambda));
er5zero=max(er5,erzero);
er5_max=max(er5zero);
er5_ave=10^(sum(log10(er5zero))/double(n));
disp(sprintf('Max: %.1e',er5_max));
disp(sprintf('Ave: %.1e',er5_ave));

disp('(5*)誤差が同じ部分と違う部分を格納');

% 異なる誤差の配列を作成
% MATLABとFORTRAN
% 誤差が同じ部分を１で違う部分を０で出力
eq = abs(er3zero-er4zero)== 0;
% 誤差が違う部分のみを出力
er3zero_not_eq = er3zero(~eq);
er4zero_not_eq = er4zero(~eq);
n1 = numel(er3zero_not_eq);
disp(sprintf('MATLABとFORTRANの誤差が異なる数: %d 個',n1));
% 誤差が同じ部分の配列のみ出力
er_zero_eq = er3zero(eq);
disp(sprintf('MATLABとFORTRANの誤差が同じ数: %d 個',numel(er_zero_eq)));

% CとFORTRAN
% 誤差が同じ部分を１で違う部分を０で出力
eq_1 = abs(er5zero-er4zero)== 0;
% 誤差が違う部分のみを出力
er5zero_not_eq = er5zero(~eq_1);
er4zero_not_eq_1 = er4zero(~eq_1);
n2 = numel(er5zero_not_eq);
disp(sprintf('CとFORTRANの誤差が異なる数: %d 個',n2));
% 誤差が同じ部分の配列のみ出力
er_zero_eq_1 = er5zero(eq_1);
disp(sprintf('CとFORTRANの誤差が同じ数: %d 個',numel(er_zero_eq)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5*)平均と分散
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('(5*)分散と平均');

% 異なる誤差の配列の平均と分散
var3 = var(er3zero_not_eq);
var4 = var(er4zero_not_eq);
disp(sprintf('MATLABの分散: %.1e',var3));
disp(sprintf('FORTRANの分散: %.1e',var4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (5**)得点計算
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -14.5を真ん中としてそれより小さいならプラス、大きいならマイナスとして計算
disp('(5**)得点を表示');
for i=1:n2
    p3(i) = -14.5-log10(er3zero_not_eq(i));
    p4(i) = -14.5-log10(er4zero_not_eq(i));
end
disp(sprintf('MATLABの得点: %g点',sum(p3)));
disp(sprintf('FORTRANの得点: %g点',sum(p4)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (6)グラフ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name_fig = append('matrix',string(x),'_','test',string(y),'_','fig');
mkdir(folder_name_fig);

% グラフ
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
%%%%% 対数の名前どうする
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

% 異なる誤差のみを出力したグラフ(MATLABとFORTRAN)
% タイトル　何のヒストグラム？　名前決める　横軸loge k?

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
histogram(log10(er5zero_not_eq));title('dqds (C):誤差が異なる');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'not_eq_C_er.fig')
saveas(gcf,'not_eq_C_er.png')
movefile('not_eq_C_er.fig',(folder_name_fig));
movefile('not_eq_C_er.png',(folder_name_fig));

figure(7);
histogram(log10(er3zero_not_eq));title('dqds (MATLAB):誤差が異なる');axis([-17 -13 0 35]);grid on;
xlabel('Relative errors ()');
ylabel('Number');
saveas(gcf,'not_eq_MATLAB_er.fig')
saveas(gcf,'not_eq_MATLAB_er.png')
movefile('not_eq_MATLAB_er.fig',(folder_name_fig));
movefile('not_eq_MATLAB_er.png',(folder_name_fig));

figure(8);
histogram(log10(er4zero_not_eq));title('dqteqr (FORTRAN):誤差が異なる');axis([-17 -13 0 35]);grid on;
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
xlabel('反復回数');
ylabel('Relative errors');
legend('dqds( C )');
saveas(gcf,'loop_er.fig')
saveas(gcf,'loop_er.png')
movefile('loop_er.fig',(folder_name_fig));
movefile('loop_er.png',(folder_name_fig));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (7)ファイルへの出力
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_name = append('matrix',string(x),'_','test',string(y));
mkdir(folder_name);
% 比較結果（実行時間、総反復回数、一個あたりの、最大値、平均、分散、同じ数と違う数）
file_name = append('data',string(y),'_','matrix','_',string(x),'.txt');
fid = fopen((file_name),'w');
% 使用したテスト行列の情報
fprintf(fid,'(1)情報\r\n');
name = append('テスト行列',string(y));
fprintf(fid,'テスト行列名:%s\n',name);
size = string(x);
fprintf(fid,'行列サイズ:%s\n',size);
fprintf(fid,'=============================================\r\n');
% 実行時間
fprintf(fid,'(2)実行時間\r\n');
fprintf(fid,'FORTRAN_Computing_time %g [s]\r\n',time_dpteqr);
fprintf(fid,'MATLAB_Computing_time %g [s]\r\n',time_matlab_dqds);
fprintf(fid,'=============================================\r\n');
% 総反復回数
fprintf(fid,'(3)総反復回数\r\n');
fprintf(fid,'MATLAB: %d\r\n',nn3);
fprintf(fid,'=============================================\r\n');
% 誤差評価
fprintf(fid,'(4)誤差評価\r\n');
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
% 誤差が同じ部分と違う部分
fprintf(fid,'(5)誤差が同じ個数と違う個数\r\n');
fprintf(fid,'誤差が異なる数: %d 個\r\n',n1);
fprintf(fid,'誤差が同じ数: %d 個\r\n',numel(er_zero_eq));
fprintf(fid,'=============================================\r\n');
% 分散
fprintf(fid,'(6)誤差が異なる部分の分散\r\n');
fprintf(fid,'MATLABの分散: %.1e\r\n',var3);
fprintf(fid,'FORTRANの分散: %.1e\r\n',var4);
fprintf(fid,'=============================================\r\n');
% 得点
fprintf(fid,'(7)得点\r\n');
fprintf(fid,'MATLABの得点: %g点\r\n',sum(p3));
fprintf(fid,'FORTRANの得点: %g点\r\n',sum(p4));
fprintf(fid,'二つの得点差: %g点\r\n',abs(sum(p3)-sum(p4)));
fprintf(fid,'=============================================\r\n');

fclose(fid);

% 真値
file_name1 = 'lambda.txt';
fid = fopen((file_name1),'w');
fprintf(fid,'固有値(真値)\r\n');
fprintf(fid,'%g\n',lambda);
fclose(fid);
% FORTRANの固有値
file_name2 = 'FORTRAN_eig.txt';
fid = fopen((file_name2),'w');
fprintf(fid,'固有値(FORTRAN)\r\n');
fprintf(fid,'%g\n',eig_dpteqr);
fclose(fid);
% MATLABの固有値
file_name3 = 'MATLAB_eig.txt';
fid = fopen((file_name3),'w');
fprintf(fid,'固有値(MATLAB)\r\n');
fprintf(fid,'%g\n',lambda3);
fclose(fid);
% Cの固有値
file_name4 = 'C_eig.txt';
fid = fopen((file_name4),'w');
fprintf(fid,'固有値(c)\r\n');
fprintf(fid,'%g\n',dqds_c);
fclose(fid);

% 誤差MATLAB
file_name5 = 'MATLAB_er.txt';
fid = fopen((file_name5),'w');
fprintf(fid,'誤差(c)\r\n');
fprintf(fid,'%g\n',er3);
fclose(fid);
% 誤差FORTRAN
file_name6 = 'FORTRAN_er.txt';
fid = fopen((file_name6),'w');
fprintf(fid,'誤差(FORTRAN)\r\n');
fprintf(fid,'%g\n',er4);
fclose(fid);
% 誤差C
file_name7 = 'C_er.txt';
fid = fopen((file_name7),'w');
fprintf(fid,'誤差(c)\r\n');
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
