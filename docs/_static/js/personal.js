    function crypto(sub) {
        const e = sub.parentNode.firstElementChild      //密文部分
        var key = sub.previousElementSibling.value  //用户输入的密钥
        while (key.length < 16){
            key = key + "\0"
        }

        const origin = e.getAttribute("data-origin")    //试探明文
        const now = e.getAttribute("data-now")          //试探密文
        const res = decrypt(now, key)
        if (res == origin) {    //如果密钥正确
            const t = e.innerHTML.replace(/\s/g, "")
            sub.parentNode.innerHTML = decrypt(t, key)
        } else {
            alert("密码错误！")
        }

    }
    // n位随机数生成
    function randomNum(n) {
        let sString = "";
        let strings = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
        for (i = 0; i < n; i++) {
            ind = Math.floor(Math.random() * strings.length);
            sString += strings.charAt(ind);
        }
        return sString
    }

    //AES-128-CBC-ZeroPaddingjie密
    function decrypt(data, key) {
        iv = CryptoJS.MD5(key).toString().substring(0, 16);  //取转化为md5格式的前面16位字符
        key = CryptoJS.enc.Utf8.parse(key);  //解析后的key
        iv = CryptoJS.enc.Utf8.parse(iv); //解析后的iv
        encrypted = CryptoJS.AES.decrypt(data, key, { //解密
            iv: iv,
            mode: CryptoJS.mode.CBC,
            padding: CryptoJS.pad.Pkcs7
        });
        return encrypted.toString(CryptoJS.enc.Utf8)

    }
    function encrypt(data, key) {
        while (key.length < 16){
            key = key + "\0"
        }

        iv = CryptoJS.MD5(key).toString().substring(0, 16);  //取转化为md5格式的前面16位字符
        key = CryptoJS.enc.Utf8.parse(key);  //解析后的key
        iv = CryptoJS.enc.Utf8.parse(iv); //解析后的iv
        encrypted = CryptoJS.AES.encrypt(data, key, { //j加密
            iv: iv,
            mode: CryptoJS.mode.CBC,
            padding: CryptoJS.pad.Pkcs7
        });
        return encrypted.toString()

}