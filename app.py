from flask import Flask, render_template, request, redirect, url_for, send_file, session
import os
import pandas as pd
import glob2 as glob
import matplotlib.pyplot as plt
import pybase64
from io import BytesIO
import shutil
import secrets

app = Flask(__name__)
# ユーザー情報
users = {
    'Proteome': {'username': 'Proteome', 'password': 'Proteome'},
}

app.config['SECRET_KEY'] = secrets.token_hex(16)

@app.route('/')
def main_page():
    return render_template('login.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        # 画面で入力された情報を取得
        username = request.form['username']
        password = request.form['password']

        # ログイン可否を判定
        if username in users and users[username]['password'] == password:
            session['username'] = username

            # ログイン成功でapp.htmlを返す
            return redirect(url_for('MSapp'))
        else:
            return render_template('login.html', error='Invalid credentials')

    # GETの場合はログイン画面へ戻す
    return render_template('login.html')


@app.route('/logout', methods=['GET', 'POST'])
def logout():
    session.pop('username', None)
    return redirect(url_for('login'))

# 操作画面
@app.route('/MSapp')
def MSapp():
    if 'username' in session:
        xlsx_files = [file for file in os.listdir('files') if file.endswith('.xlsx')]
        return render_template('MSapp.html', xlsx_files=xlsx_files)
    else:
        return redirect(url_for('login'))

# ファイルのアップロードを行う
@app.route('/upload', methods=['POST'])
def upload():
    files = request.files.getlist('file')
    for i, file in enumerate(files):
        file_name = file.filename
        file_path = os.path.join('files', file_name)
        file.save(file_path)
    return redirect(url_for('MSapp'))

# 平均したファイルのダウンロードを行う
@app.route('/download/')
def download():
    return send_file('result/YYMMDD_sample_XX-XX_ave.xlsx', mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     as_attachment=True)

#グラフの画像のダウンロードを行う
@app.route('/download_img/')
def download_img():
    return send_file('result/YYMMDD_sample_XX-XX.png', as_attachment=True)

# ペプチドタグの全データファイルのダウンロードを行う
@app.route('/download_ratio/')
def download_ratio():
    return send_file('result/YYMMDD_sample_XX-XX_ratio_all.xlsx', mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     as_attachment=True)

# ペプチドタグの最大値のみのファイルのダウンロードを行う
@app.route('/download_ratio_max/')
def download_ratio_max():
    return send_file('result/YYMMDD_sample_XX-XX_ratio_max.xlsx', mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                     as_attachment=True)

# ファイルの削除を行う
@app.route('/delete/<string:file>')
def delete(file):
    delete_file_path = os.path.join('files', file)
    os.remove(delete_file_path)
    return redirect(url_for('MSapp'))

#ファイルの一括削除
@app.route('/delete_all/')
def del_all():
    shutil.rmtree('files/')
    os.mkdir('files')
    return redirect(url_for('MSapp'))


# グラフ作成を行う
@app.route('/graph/', methods=['POST'])
def make_graph():
    path = "files"
    rfiles = glob.glob(path + "*.xlsx")
    title = request.form['title']
    ver = request.form['ver']
    intens = request.form['intens']
    remove_ker = request.form['remove_ker']

    rfile_list = []
    for filename in rfiles:
        if ver == "PD13":
            df = pd.read_excel(filename, 'Pre_Q', index_col=None)
            df = df[['GeneName', 'Sequence', 'Modifications', 'Charge', 'Precursor Area', 'RT [min]']]
            rfile_list.append(df)
        elif ver == "PD22":
            df = pd.read_excel(filename, 'Pre_Q', index_col=None)
            df = df[['GeneName', 'Annotated Sequence', 'Modifications', 'Charge', 'Precursor Abundance', 'RT [min]']]
            rfile_list.append(df)

    df = pd.concat(rfile_list, axis=0, ignore_index=True)
    df = df.reset_index(drop=True)
    df = df.fillna('NaN')

    if remove_ker == 'yes':
        ker_index = df.index[df['GeneName'].str.match('\(.*')]
        df = df.drop(ker_index)

    if ver == "PD13":
        df_ave = df.groupby(['GeneName', 'Sequence', 'Modifications', 'Charge'], as_index=False).mean()
    elif ver == "PD22":
        df_ave = df.groupby(['GeneName', 'Annotated Sequence', 'Modifications', 'Charge'], as_index=False).mean()


    df_ave.to_excel("result/YYMMDD_sample_XX-XX_ave.xlsx")

    tag_df = df_ave[df_ave['GeneName'].str.match("^\w+_PepTag$", na=False)]

    x = df_ave.loc[:, 'RT [min]']
    if ver == "PD13":
        y = (df_ave.loc[:, 'Precursor Area']) / int(intens)
    elif ver == "PD22":
        y = (df_ave.loc[:, 'Precursor Abundance']) / int(intens)
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.scatter(x, y, facecolor='None', edgecolors='grey', marker='o', label='other peptide', alpha=0.5)

    i = 0
    if ver == "PD13":
        for item in zip(tag_df['GeneName'], tag_df['Modifications'], tag_df['Charge'],
                        tag_df['Precursor Area'], tag_df['RT [min]']):
            if (item[1] == 'NaN') & (item[2] == 2):
                ax.scatter(item[4], item[3] / 10000, c=[plt.cm.tab20(i)], s=120, marker='D',
                       label='tag(' + item[0][:4] + ')')
                i = i + 1
            elif item[2] == 2:
                ax.scatter(item[4], item[3] / 10000, c=[plt.cm.tab20(i)], s=120, marker='D',
                       label='tag(' + item[0][:4] + ')_heavy')
                i = i + 1
    elif ver == "PD22":
        for item in zip(tag_df['GeneName'], tag_df['Modifications'], tag_df['Charge'],
                        tag_df['Precursor Abundance'], tag_df['RT [min]']):
            if (item[1] == 'NaN') & (item[2] == 2):
                ax.scatter(item[4], item[3] / 10000, c=[plt.cm.tab20(i)], s=120, marker='D',
                        label='tag(' + item[0][:4] + ')')
                i = i + 1
            elif item[2] == 2:
                ax.scatter(item[4], item[3] / 10000, c=[plt.cm.tab20(i)], s=120, marker='D',
                        label='tag(' + item[0][:4] + ')_heavy')

    ax.set_yscale('log')
    ax.set_xlim(0, 160)
    ax.set_ylim(1, 10 ** 7)
    ax.set_title(title, fontsize=40)
    ax.set_xlabel('RT (min)', fontsize=32)
    ax.set_ylabel('Intensity(/' + intens + ')', fontsize=32)
    ax.tick_params(labelsize=30)
    ax.set_axisbelow(True)
    ax.grid(c='gray')
    plt.legend(bbox_to_anchor=(1.5, 1), loc='upper right', borderaxespad=0.5, fontsize=20)
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    plt.tight_layout()
    fig.savefig("result/YYMMDD_sample_XX-XX.png")

    buf = BytesIO()
    fig.savefig(buf, format='png')

    data = pybase64.b64encode(buf.getbuffer()).decode('ascii')

    return render_template(
        'graph.html', data=data, table=(tag_df.to_html(classes='mystyle'))
    )

@app.route('/ratio/')
def ratio():
    path = "files"
    rfiles = glob.glob(path + "*.xlsx")

    file_list = []
    for filename in rfiles:
        df = pd.read_excel(filename, 'Pre_Q', index_col=None)
        df = df[['GeneName', 'Sequence', 'Modifications', 'Charge', 'Precursor Area', 'RT [min]', 'Light/Heavy', 'Light', 'Heavy', 'Spectrum File']]
        df = df[df['GeneName'].str.match("^\w+_PepTag$", na=False)]
        file_list.append(df)

    df = pd.concat(file_list, axis=0, ignore_index=True)
    df = df.reset_index(drop=True)
    df = df.fillna('NaN')
    df = df.sort_values(['GeneName', 'Precursor Area'], ascending=[True, False])
    df = df.reset_index(drop=True)
    tag_df = df[df['Charge'] == 2]
    tag_df = tag_df.drop_duplicates(subset='GeneName', keep='first')

    df.to_excel("result/YYMMDD_sample_XX-XX_ratio_all.xlsx")
    tag_df.to_excel("result/YYMMDD_sample_XX-XX_ratio_max.xlsx")
    return render_template(
        'ratio.html', table=(df.to_html(classes='mystyle'))
    )


# readme
@app.route('/readme')
def readme():
    return render_template('readme.html')

if __name__ == '__main__':
    app.run(debug=True)



# エラーメッセージ
@app.errorhandler(500)
def internal_error(error):
    return render_template('500.html', msg="入力に誤りがあります"), 500