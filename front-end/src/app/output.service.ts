import { BehaviorSubject, Observable } from 'rxjs';
import { DomSanitizer } from '@angular/platform-browser';
import { Injectable } from '@angular/core';
import { MatSnackBar } from '@angular/material/snack-bar';

import { BackendSocketClient } from './backend-socket.client';
import { Block } from './block.interface';
import { NewBlock, Request } from './request';
import { ImageInfo, Output } from './output';
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

  constructor(private userIdService: UserIdService, private backendSocketClient: BackendSocketClient, private snackBar: MatSnackBar, private sanitizer: DomSanitizer) {

    this.userIdService.userId.subscribe(
      (res) => { this.userId = res; },
    );

    // TODO should this just be a normal http request? do we really need socket io lib?
    // e.g. https://angulargems.beehiiv.com/p/web-sockets-in-angular
    this.backendSocketClient.listen((msg: string) => {
      const res = JSON.parse(msg);
      if (res.text) {
        const processedResponse: Output = {
          blockId: res.blockId,
          text: res.text
        };
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } else if (res.image) {
        const imageString = res.image.image as string;
        const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
        const objectURL = 'data:image/png;base64,' + processedString;
        const newImage = this.sanitizer.bypassSecurityTrustUrl(objectURL);
        const processedResponse: Output = {
          blockId: res.blockId,
          image: {
            image: newImage,
            altText: res.image.alt_text
          }
        };
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } else if (res.image_list) {
        const imageList: ImageInfo[] = [];
        for (const imageListItem of res.image_list) {
          const imageString = imageListItem.image as string;
          const processedString = imageString.substring(2, imageString.length-3).replace(/\\n/g, '');
          const objectURL = 'data:image/png;base64,' + processedString;
          const newImage = this.sanitizer.bypassSecurityTrustUrl(objectURL);
          imageList.push({
            image: newImage,
            altText: imageListItem.alt_text
          });
        }
        const processedResponse: Output = {
          blockId: res.blockId,
          imageList: imageList
        };
        const outputs = this.outputs$.getValue();
        outputs.push(processedResponse);
        this.outputs$.next(outputs);
      } else if (res.end_connection) {
        this.executingBlocks$.next(false);
      } else if (res.error) {
        this.snackBar.open(res.error, 'Close', { duration: 5000 });
        this.executingBlocks$.next(false);
      } else {
        this.snackBar.open('Received a bad response. Please refresh the page and try again.', 'Close', { duration: 5000 });
        this.executingBlocks$.next(false);
      }
    });
  }

  resetOutputs(): void {
    this.outputs$.next([]);
  }

  executeBlocks(blocks: Block[]): void {
    if (this.userId === null) {
      this.snackBar.open('No User ID. Please refresh the page and try again.', 'Close', { duration: 5000 });
    } else {
      this.executingBlocks$.next(true);
      const newBlocksArray: NewBlock[] = [];
      for (let i = 0; i < blocks.length; i++) {
        const block: NewBlock = {};
        block['block_id'] = blocks[i].blockId;
        for (let j = 0; j < blocks[i].parameters.length; j++) {
          block[blocks[i].parameters[j].key] = blocks[i].parameters[j].value;
        }
        newBlocksArray.push(block);
      }
      const request: Request = {
        'user_id': this.userId,
        'blocks': newBlocksArray,
      };
      this.backendSocketClient.sendRequest(request);
    }
  }
}
