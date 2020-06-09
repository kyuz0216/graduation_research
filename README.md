# graduation_research
graduation research B4  

# gitメモ  

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
