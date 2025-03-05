#!/usr/bin/env python3

import subprocess
import os

class FileManager:
    def __init__(self, remote:str, local_master:str, cloud_master:str):
        """Initialize attributes"""
        self.remote = remote
        self.local_master = local_master
        self.cloud_master = cloud_master

    def convertCloudToLocal(self, filename:str) -> str:
        """converts cloud location into local path"""
        cloud = f'{self.remote}:/{self.cloud_master}/'
        file_path = filename.replace(cloud, '')
        local = f'{self.local_master}/{file_path}'
        return local

    def convertLocalToCloud(self, filename:str) -> str:
        """converts local location into cloud location"""
        local = f'{self.local_master}/'
        file_path = filename.replace(local, '')
        cloud = f'{self.remote}:/{self.cloud_master}/{file_path}'
        return cloud

    def uploadData(self, filename:str):
        """uploads local file to cloud"""

        cloud_loc = self.convertLocalToCloud(filename)
        upload = subprocess.run(['rclone', 'copy', filename, cloud_loc], capture_output = True, text = True)
        if upload.returncode == 0:
            print(f"Upload successful! File uploaded to {cloud_loc}")
        else:
            print(f"Upload failed. Error:\n{upload.stderr}")

    def downloadData(self, filename:str):
        """downloads cloud file to local path"""

        local_loc = self.convertCloudToLocal(filename)
        download = subprocess.run(['rclone', 'copyto', filename, local_loc], capture_output = True, text = True)
        if download.returncode == 0:
            print(f"Download successful! File downloaded to {local_loc}")
        else:
            print(f"Download failed. Error:\n{download.stderr}")

    def __str__(self) -> str:
        return f'Remote: {self.remote} \nLocal master path: {self.local_master} \nCloud master path: {self.cloud_master}'


fm = FileManager("dropbox_remote", "/home/anagha/BIOL8802/topic4/hw2", "Anagha Mohana Krishna/dropbox_rclone")
print(fm)
fm.uploadData('/home/anagha/BIOL8802/topic4/hw2/local_rclone/test.txt')
fm.downloadData('dropbox_remote:/Anagha Mohana Krishna/dropbox_rclone/test1/Biopython_re_comparision.png')