import { BehaviorSubject, Observable } from 'rxjs';
import { DomSanitizer } from '@angular/platform-browser';
import { Injectable, isDevMode } from '@angular/core';
import { Socket, io } from 'socket.io-client';
import { MatSnackBar } from '@angular/material/snack-bar';

import { Block } from './block.interface';
import { Output } from './output';
import { OutputServiceInterface } from './output.service.interface';
import { UserIdService } from './user-id.service';

@Injectable({
  providedIn: 'root',
})
export class OutputService implements OutputServiceInterface {
  private readonly executingBlocks$: BehaviorSubject<boolean> = new BehaviorSubject<boolean> (false);
  readonly executingBlocks: Observable<boolean> = this.executingBlocks$.asObservable();

  private readonly outputs$: BehaviorSubject<Output[]> = new BehaviorSubject<Output[]> ([]);
  readonly outputs: Observable<Output[]> = this.outputs$.asObservable();

  private userId: string | null = null;

  private socket: Socket;

  constructor(private userIdService: UserIdService, private snackBar: MatSnackBar, private sanitizer: DomSanitizer) {
    this.userIdService.userId.subscribe(
      (res) => { this.userId = res; },
    );
    const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
    // TODO: Are these timeouts required?
    this.socket = io(url, {
      autoConnect: false,
      timeout: 1000000,   // 1000 seconds,
      ackTimeout: 1000000 // 1000 seconds,
    });

    this.socket.on('json', (msg: string) => {
      const res = JSON.parse(msg);
      console.log('res', typeof msg, msg);
      if (res.text) {
        const processedResponse: Output = {};
        processedResponse.text = res.text;
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } else if (res.img && res.alttext) {
        const imageString = res.img as string;
        const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
        const objectURL = 'data:image/png;base64,' + processedString;
        const newImg = this.sanitizer.bypassSecurityTrustUrl(objectURL);
        const processedResponse: Output = {};
        processedResponse.img = newImg;
        processedResponse.alttext = res.alttext;
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } else if (res.end_connection) {
        this.executingBlocks$.next(false);
      } else if (res.error) {
        this.snackBar.open(res.error, 'Close', { duration: 5000 });
        this.executingBlocks$.next(false);
      } else {
        this.snackBar.open('Received a bad response, please refresh the page and try again', 'Close', { duration: 5000 });
        this.executingBlocks$.next(false);
      }
    });

    this.socket.connect();
  }

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  executeBlocks(blocks: Block[]): void {
    if (this.userId === null) {
      this.snackBar.open('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    } else {
      this.executingBlocks$.next(true);
      type NewBlock = { [key: string]: string | number };
      const newBlocksArray: NewBlock[] = [];
      for (let i = 0; i < blocks.length; i++) {
        const block: NewBlock = {};
        block['block_id'] = blocks[i].blockId;
        for (let j = 0; j < blocks[i].parameters.length; j++) {
          block[blocks[i].parameters[j].key] = blocks[i].parameters[j].value;
        }
        newBlocksArray.push(block);
      }    
      type Message = { [key: string]: string | NewBlock[] };
      const message: Message = {};
      message['user_id'] = this.userId;
      message['blocks'] = newBlocksArray;
      this.socket.emit('json', message);
    }
  }
}
