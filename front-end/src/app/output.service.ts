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
    let url = '';
    if (isDevMode()) {
      url = 'ws://127.0.0.1:5000';
    } else {
      url = '';
    }
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

      const blockStringArray: string[] = [];
      for (let i = 0; i < blocks.length; i++) {
        const block: Block = blocks[i];
        let blockString: string = block.blockId;
        for (let j = 0; j < block.parameters.length; j++) {
          blockString = blockString.concat(',');
          blockString = blockString.concat(block.parameters[j].key);
          blockString = blockString.concat('=');
          blockString = blockString.concat(block.parameters[j].value.toString());
        }
        blockStringArray.push(blockString);
      }
      const blockStrings: ReadonlyArray<string> = blockStringArray;
      type Params = { [key: string]: string | ReadonlyArray<string>};
      const params: Params = {};
      params['user_id'] = this.userId;
      params['block_strings'] = blockStrings;

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
          }
          const outputs = this.outputs$.getValue();
          outputs.push(processedResponse);
          this.outputs$.next(outputs);
        }
      );

      this.subject.next(params);

      // try {
      //   null;
      // } catch (e) {
      //   if (e instanceof HttpErrorResponse && e.status == 406) {
      //     this.snackBar.open('The blocks you have executed are not a valid order. Please check the blocks and try again.', 'Close', { duration: 5000 });
      //   } else if (e instanceof HttpErrorResponse && e.status == 400 && e.error) {
      //     this.snackBar.open(e.error.description, 'Close', { duration: 5000 });
      //   } else {
      //     this.snackBar.open('There has been an error. Please refresh the page and try again.', 'Close', { duration: 5000 });
      //   }
      // }

    }
  }
}
