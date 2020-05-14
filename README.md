# LAMMPS_nonlinear_acoustic_wave

outerview
check num of arrage
if you do not need visualization, delete 'dump'
Use OVITO to check molphorogy
FFT configuration check
FFT is automized using a python file
arraving time should be defined consistnatly
change x_detec, A2, Time to get beta

makebcclat.py
>FeまたはFe-Cuを並べる。
(R ,type) type= vacancy:1 , Cu: 2
and RATE or num (=N_X,N_Y)are all parameters
python makebcc.py
とコマンド入力

Fe.lat
>makebcclat.pyの出力ファイル。粒子数などを確認

in.lammps
>lammpsのinputファイル
系によって
pair_coeff  * * ./new_potential/FeCuNi.eam.alloy
dump
の末尾を変えること
Fe系のとき粒子グループは4つ 
Fe Fe Fe Fe 
Fe-Cu系のとき粒子グループは5つ
Fe Fe Fe Fe Cu
可視化したいときはdumpをつける（2行）

FeCuNialloy.eam
>potentialファイル

move.sh
>スパコンを動かすファイル
sbatch move.sh
とコマンド入力

log.lammps
>lammpsの出力ファイル。必要な出力はこっちではない

outp.txt
>こちらが私が設定した出力
time source_loc detec_loc temperature volume energy の順
source_loc 左端の原子群の平均位置
detec_loc　波動の伝搬を読み取るために中央当たりの原子群を指定。その粒子群の平均位置

スペースを,に変換した後
python読み込み用のcsvファイルを作る
cat unti.txt | tr " " "," > outp.csv 
とコマンド入力
FFT_oritin.py
> python FFT_origin.py
とコマンド入力すると必要な情報と非線形パラメータβが出てくる。
非線形超音波の可視化までできるが、スパコンでは可視化はできないのでファイルをスパコンの外に出して実行してみると、波がきちんと可視化されていることを確認できる。

