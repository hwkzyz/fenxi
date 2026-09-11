%使用该代码前根据所使用的采集卡参数进行修改代码，需要修改①文件存储路径；②头文件长度；③采集卡位数；④量程；⑤通道数量。
clc;clear;
fip=fopen('D:\12345.bin');%①文件存储路径；
fseek(fip,0,'eof');
fsize=ftell(fip); 
datacount=(fsize-0)/2;%头文件长度为0（ART-SCOPE面板保存的数据头文件长度为0）%②头文件长度；
startPosition = 0; 


fseek(fip, startPosition, 'bof');% 读取指定的字节数
c1=fread(fip,inf,'uint16');
rows=bin2dec('101');
rows=size(c1,1);
C2=ones(rows,1)*4095;% ③采集卡位数12位 0xFFF（4095），14位0x3FFF（16383），16位0xFFFF（65535）
c3=bitand(c1,C2);       
fclose(fip);

data1 = (10000/4096).*c3-5000; 
%④(10000/4096).*c3-5000该公式表示使用12位（4096），±5V量程采集卡保存到的数据由码值转成实际电压值（mV）；10000对应±5V（满量程10000mV）（±1V：2000），4096对应12位采集卡（14bit 16384，16bit 65536），5000对应5V；

data2=reshape(data1,[1,rows/1]);     %⑤[1,rows/1]指采集通道为1，[2,rows/2]指采集通道为2，以此类推，根据实际设置；


% 获取矩阵的行数
numRows = size(data2, 1);%1是按行绘图，2是按列绘图

% 创建一个图形窗口
figure;

% 循环遍历每一行数据并绘图
for i = 1:numRows
% 提取当前行的数据
rowData = data2(i,:);%按行这样写，按列调换括号内 内容的位置

% 绘制当前行的数据
plot(rowData, 'DisplayName', ['Row ', num2str(i)]);

% 保持当前图形，以便在同一图中绘制多行数据
hold on;
end

% 设置图形属性
title('Plot of Each Row in the Matrix');
xlabel('Data Point Index');
ylabel('Value');
legend;

% 自动调整坐标轴范围
axis tight;

% 关闭图形保持状态
hold off;





