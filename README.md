# graduation_research
graduation research B4  

# 全体の流れ
## ファイルを変更した時
1.ローカル(MATLAB)でコードを編集する。  
2.ブランチを作成する  
3.正しく動いていれば、その変更内容をresearchフォルダのmanin_new.mに反映する。（めんどくさいけど、コピペか置き換えか。）  
  ＊この際、manin_new.m以外に追加のファイルがあれば、それもresearchフォルダに入れる。  
    その場合、その有無をPull Requestの文のところに書いておく。  
4.git add→git commit ~　→git pushでgithubに反映する  
5.githubにブラウザでアクセスし、Pull Requestを出す。  対面ではない時は、詳しくPull Requestの概要欄に詳細を記入する  
6.相手にPull Requestの承認(マージ)をしてもらう。(自分でもボタンは押せるが、相手にしてもらったほうが良い。)  

## 相手の変更を受け取りたい時
1.masterブランチにする  
2.git pullをする  
3.変更部分を自分のローカル(matlab)のファイルに反映させる  

### ブランチのルール
ブランチ名は'自分の名前_番号'にする。  
### ブランチの作り方
```
git checkout -b ブランチ名
```
### 現在のブランチの確認
```
git branch
```

### ブランチの切り替え
```
git checkout ブランチ名
```

### Pull Requestを送るまで
```
git add .
git commit -m "変更内容の説明"
git push origin ブランチ名
```
addできているかの確認  
```
git status
```

### PULLする方法
masterブランチに切り替えてから
```
git pull origin master
```
