import os
import sys

from bs4 import BeautifulSoup
import base64
import hashlib
from Crypto.Util.Padding import pad
from Crypto.Cipher import AES

key = os.environ["MYKEY"]

while(len(key.encode()) <16):
    key += '\0'

class EnAesCrypt:
    """
    AES-128-CBC解密
    """

    def __init__(self, data, key):
        """
        :param data: 加密后的字符串
        :param key: 随机的16位字符
        :param pad: 填充方式
        """
        self.key = key.encode("utf-8")
        self.data = data.encode("utf-8")

        hash_obj = hashlib.md5()  # 构造md5对象
        hash_obj.update(self.key)  # 进行md5加密,md5只能对byte类型进行加密
        res_md5 = hash_obj.hexdigest()  # 获取加密后的字符串数据
        self.iv = res_md5[:16].encode("utf-8")

    def encrypt_aes(self):
        """AES-128-CBC加密"""

        my_aes = AES.new(self.key, AES.MODE_CBC, self.iv)
        ct_bytes = my_aes.encrypt(pad(self.data, AES.block_size))
        ct = base64.b64encode(ct_bytes).decode('utf-8')
        return ct


def encode_html(file):

    soup = BeautifulSoup(open(file,"r"), 'html.parser')
    encode_soup = soup.find_all(id="main-content")
    encode_tag = encode_soup[0]

    data_string = ''
    for child in encode_tag.children:
        try:
            data_string += child.prettify()
        except:
            data_string += str(child)


    encode_tag.clear()

    # encode
    encode = EnAesCrypt(data_string, key)
    encode_data_str = encode.encrypt_aes()

    data_origin = "stDkxwE"
    data_now = EnAesCrypt(data_origin, key).encrypt_aes()

    append_str = f'''
    <div>
    <div data-origin="{data_origin}" data-now="{data_now}" style="display: none" class="cryptoText">
        {encode_data_str}
    </div>
    请输入密码查看隐藏内容：<input type="password" />
    <input type="submit" onclick="crypto(this)" /> 
    </div>
    '''

    encode_tag.append(BeautifulSoup(append_str, 'html.parser').div)

    with open(file, 'w') as w:
        w.write(soup.prettify())

if __name__ == "__main__":
    encode_html(sys.argv[1])
