import { BehaviorSubject, Observable } from 'rxjs';
import { DomSanitizer } from '@angular/platform-browser';
import { Injectable, isDevMode } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';
import { WebSocketSubject, webSocket } from 'rxjs/webSocket';

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

  private subject: WebSocketSubject<unknown>;

  constructor(private userIdService: UserIdService, private snackBar: MatSnackBar, private sanitizer: DomSanitizer) {
    this.userIdService.userId.subscribe(
      (res) => { this.userId = res; },
    );
    const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
    this.subject = webSocket(url);
  }

  resetOutputs(): void {
    this.outputs$.next([]);
  }
  
  executeBlocks(blocks: Block[]): void {
    if (this.userId === null) {
      this.snackBar.open('No User ID, please refresh the page and try again', 'Close', { duration: 5000 });
    } else {
      this.executingBlocks$.next(true);

      this.subject.subscribe(
        (msg: any) => {
          const processedResponse: Output = {};
          if (msg.text) {
            processedResponse.text = msg.text;
          } else if (msg.img && msg.alttext) {
            const imageString = msg.img as string;
            const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
            const objectURL = 'data:image/png;base64,' + processedString;
            const newImg = this.sanitizer.bypassSecurityTrustUrl(objectURL);
            processedResponse.img = newImg;
            processedResponse.alttext = msg.alttext;
          } else if (msg.end_connection) {
            this.subject.complete();
            this.executingBlocks$.next(false);
          } else if (msg.error) {
            this.snackBar.open(msg.error, 'Close', { duration: 5000 });
            this.subject.complete();
            this.executingBlocks$.next(false);
          } else {
            this.snackBar.open('Received a bad response, please refresh the page and try again', 'Close', { duration: 5000 });
            this.subject.complete();
            this.executingBlocks$.next(false);
          }
          const outputs = this.outputs$.getValue();
          outputs.push(processedResponse);
          this.outputs$.next(outputs);
        }
      );

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
      this.subject.next(message);
    }
  }
}
