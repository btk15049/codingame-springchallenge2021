# FAQ

## 自分用の備忘録

### ssh 系がうまく動かない場合

ssh-agent を起動する

bash の場合

```bash
eval $(ssh-agent)
ssh-add
```

fish の場合

```fish
eval (ssh-agent -c)
```

### cargo 実行時に `error authenticating: no auth sock variable; class=Ssh` と出る

上記の ssh-agent の対応をする
