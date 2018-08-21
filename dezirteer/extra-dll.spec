# -*- mode: python -*-

block_cipher = None


a = Analysis(['pushkin\\AppData\\Roaming\\Python\\Python36\\site-packages\\scipy\\extra-dll', 'main_gui.py'],
             pathex=['C:\\Users\\alex', 'C:\\odrive\\Amazon Cloud Drive\\cloud\\Developing\\dezirteer\\dezirteer'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='extra-dll',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=True )
